import polars as pl
import numpy as np
from pathlib import Path
import cmdstanpy as cmd
import arviz as az
import matplotlib.pyplot as plt


def main():
    here = Path(__file__).parent.resolve()
    root = here.parent
    data_dir = root / "data/raw"
    model = cmd.CmdStanModel(stan_file= here / "fluxomics.stan")
    data_s1 = pl.read_csv(data_dir / "250218_Fluxomics_split1_900uL_g3p_excluded.csv")
    data_s2 = pl.read_csv(data_dir / "250218_Fluxomics_split2_900uL.csv")
    data = pl.concat([data_s1,data_s2])
    theoretical = pl.read_csv(data_dir / "theoretical.csv")
    undilted_samples = data.filter(pl.col("Sample Name").str.contains("_1x_"))
    theoretical = theoretical.filter(~pl.col("ID").str.contains_any(["6pgc", "oaa", "g3p"]))
    undilted_samples = undilted_samples.filter(~pl.col("Component Group Name").str.contains_any(["6pgc", "oaa", "g3p"]))
    theoretical = theoretical.with_columns(
        met = pl.col("ID").map_elements(
            lambda x: x[:-3], return_dtype=pl.String
        ),
        parent_group = pl.col("ID").map_elements(
            lambda x: x.split("_")[-1][1], return_dtype=pl.String
        ).cast(pl.Int8) + 1,
    )
    met_group = {met: ix + 1 for ix, met in enumerate(theoretical["met"].unique())}
    theoretical = theoretical.with_columns(
        met_group = pl.col("met").replace_strict(met_group),
        Theoretical = pl.col("Theoretical")/100,
    )
    undilted_samples = undilted_samples.with_columns(
        Height = pl.col("Height").cast(
            pl.Float64, strict=False
        ),
        split = pl.col("Sample Name").map_elements(
            lambda x: x.split("_")[-2], return_dtype=pl.String
        ),
        injection = pl.col("Sample Name").map_elements(
            lambda x: x.split("_")[-1], return_dtype=pl.String
        ),
        sample = pl.col("Sample Name").map_elements(
            lambda x: x.split("_")[2], return_dtype=pl.String
        ).alias("sample"),
        parent_group = pl.col("Component Name").map_elements(
            lambda x: x.split("_")[-1].split("-")[0][1], return_dtype=pl.String
        ).cast(pl.Int8) + 1
    )
    undilted_samples = undilted_samples.with_columns(
        sample_id = pl.concat_str(["split", "injection", "sample"])
    )
    undilted_samples = undilted_samples.drop_nans("Height")
    sample_group = {sample: ix + 1 for ix, sample in enumerate(undilted_samples["sample_id"].unique())}
    undilted_samples = undilted_samples.with_columns(
        met_group = pl.col("Component Group Name").replace_strict(met_group),
        sample_group = pl.col("sample_id").replace_strict(sample_group)
    )
    undilted_samples = undilted_samples.select(["Height", "parent_group", "met_group", "sample_group"])
    undilted_samples = undilted_samples.group_by(["parent_group", "met_group", "sample_group"]).sum()
    processed_samples = undilted_samples.join(theoretical, on=["parent_group", "met_group"])
    processed_samples = processed_samples.with_columns(
        met_sample = pl.concat_str(
            [pl.col("met_group"), pl.col("sample_group")],
            separator="-"
        )
    )
    sample_met_id = processed_samples["met_sample"].unique()
    sample_met_id_map = {id: ix+1 for ix, id in enumerate(sample_met_id)}
    processed_samples = processed_samples.with_columns(
        met_sample_ix = pl.col("met_sample").replace_strict(sample_met_id_map)
    )
    stan_input = {
        "N_meas": processed_samples.shape[0],
        "N_mets": int(len(met_group)),
        "N_samples": int(len(sample_group)),
        "N_met_samples": int(len(sample_met_id_map)),
        "mets": processed_samples.select(["met_group"]).to_numpy().T[0],
        "met_samples": processed_samples.select(["met_sample_ix"]).to_numpy().T[0],
        "p_theoretical": processed_samples.select(["Theoretical"]).to_numpy().T[0],
        "y": processed_samples.select("Height").to_numpy().T[0],
    }
    pf = model.pathfinder(
        data=stan_input
    )
    fit = model.sample(
        data=stan_input,
        inits=pf.create_inits(),
        chains=4,
        iter_warmup=2000,
        iter_sampling=1000,
        adapt_delta=0.99,
    )
    infd = az.from_cmdstanpy(
        posterior=fit,
        observed_data={
            "rescaled_estimated_height": processed_samples.select("Height").to_numpy().T[0]
        },
        coords={
            "met": list(met_group.keys()),
            "met_sample": list(sample_met_id_map.keys()),
        },
        dims={
            "ln_offset": ["met"],
            "true_offset": ["met"],
            "height_T": ["met_sample"],
        }
    )
    fig, ax = plt.subplots(1,1)
    az.plot_forest(infd, var_names=["ln_offset"], ax=ax)
    plt.savefig("ln_offsets.png")
    plt.close()
    fig, ax = plt.subplots(1,1)
    az.plot_forest(infd, var_names=["true_offset"], ax=ax)
    plt.savefig("offsets.png")
    plt.close()
    fig, ax = plt.subplots(1,1)
    az.plot_forest(infd, var_names=["height_T"], ax=ax)
    ax.set_xscale("log")
    plt.savefig("height_T.png")
    fig, ax = plt.subplots(1,1, figsize=(15,15))
    az.plot_forest(infd, var_names=["rescaled_height_offset"], combined=True, ax=ax)
    plt.savefig("estimated_height.png")


if __name__ == "__main__":
    main()
