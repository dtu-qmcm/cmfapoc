"""Code for preparing data."""

from collections.abc import Sequence
from pathlib import Path

from pandas.io.formats.printing import pprint_thing_encoded
from polars._typing import IntoExpr
import polars as pl

from cmfapoc.simulation import simulate

ROOT = Path(__file__).parent.parent.parent
RAW_DIR = ROOT / "data" / "raw"
OUT_DIR = ROOT / "data" / "prepared"
RAW_FILES = {
    "sergi": (RAW_DIR / "PivotTable_freqs 1.csv",),
    "daria": (
        RAW_DIR / "250218_Fluxomics_split1_900uL_g3p_excluded.csv",
        RAW_DIR / "250218_Fluxomics_split2_900uL.csv",
        RAW_DIR / "theoretical.csv",
    ),
}
SIM_ERROR_SD = 0.1


def close(expr: pl.Expr, over: IntoExpr | Sequence[IntoExpr]):
    """Close expression by sum over given columns.

    Args:
        expr: Expression to normalise.
        over: Columns to normalise over.

    Example:

    >>> df = pl.DataFrame(
        {
            "a": [1, 2, 3],
            "b": ["x", "x", "y"],
        }
    )
    >>> df.with_columns(c=close(pl.col("a"), over=["b"]))
        ┌─────┬─────┬──────────┐
        │ a   ┆ b   ┆ c        │
        │ --- ┆ --- ┆ ---      │
        │ i64 ┆ str ┆ f64      │
        ╞═════╪═════╪══════════╡
        │ 1   ┆ x   ┆ 0.333333 │
        │ 2   ┆ x   ┆ 0.666667 │
        │ 3   ┆ y   ┆ 1.0      │
        └─────┴─────┴──────────┘

    """
    return expr / expr.sum().over(over)


def prepare_data_sergi(raw: pl.DataFrame) -> pl.DataFrame:
    """Prepare data."""
    out = (
        raw.rename(
            {
                "component_name": "isotopologue",
                "component_group_name": "species",
            },
        )
        .filter(pl.col("meta_value") == "peak_apex_int")  # pyright: ignore[reportUnknownMemberType]
        .drop("meta_value")
        .unpivot(
            index=["isotopologue", "species"],
            variable_name="measurement_id",
            value_name="raw_fraction",
        )
        .with_columns(
            sample_id=pl.col("measurement_id")
            .str.extract(r"ID_(\d+)_")
            .cast(int),
            replicate_id=pl.col("measurement_id")
            .str.extract(r"ID_\d+_rep_(\d+)")
            .cast(int),
            fraction=close(
                pl.col("raw_fraction"),
                over=["measurement_id", "species"],
            ),
        )
        .with_columns(is_c12=pl.col("sample_id") == 0)
    )
    for transformation in ("alr", "clr", "ilr"):
        sim = simulate(
            dataset=out,
            transformation=transformation,
            error_sd=SIM_ERROR_SD,
        )
        out = out.with_columns(sim.alias(f"sim_fraction_{transformation}"))
    return out


def prepare_data_daria(
    raw1: pl.DataFrame,
    raw2: pl.DataFrame,
    natural: pl.DataFrame,
) -> pl.DataFrame:
    """Prepare Daria's calibration data."""
    outcols = [
        "sample",
        "metabolite",
        "isotopologue",
        "measurement",
        "measured_fraction",
        "natural_fraction",
    ]
    metabolites_to_exclude = ["6pgc", "oaa", "g3p"]

    filter_natural = ~pl.col("ID").str.contains_any(metabolites_to_exclude)
    filter_msts = pl.col("Sample Name").str.contains("_1x_") & ~pl.col(
        "Component Group Name"
    ).str.contains_any(metabolites_to_exclude)
    natural = (
        natural.filter(filter_natural)
        .with_columns(
            isotopologue=pl.col("ID").str.extract(".*_(m\\d+)$"),
            metabolite=pl.col("ID").str.extract("(.*)_m\\d+$"),
        )
        .with_columns(
            natural_fraction=close(pl.col("Theoretical"), over="metabolite")
        )
    )
    raw = pl.concat([raw1, raw2])
    new_names = {
        "Height": "measurement",
        "Sample Name": "sample",
        "Component Name": "component",
        "Component Group Name": "metabolite",
        "Sample Type": "sample_type",
    }
    msts = (
        raw.filter(filter_msts)
        .rename(new_names)
        .with_columns(
            measurement=pl.col("measurement").cast(pl.Float64, strict=False),
            split=pl.col("sample").str.extract("split(\\d)"),
            injection=pl.col("sample").str.extract("inj(\\d)"),
            qc_number=pl.col("sample").str.extract("(QC\\d?)"),
            isotopologue=pl.col("component").str.extract("(m\\d)"),
        )
        .group_by("sample", "metabolite", "isotopologue")
        .agg(pl.col("measurement").sum())
    )
    natural_with_samples = natural.join(
        msts[["metabolite", "sample"]].unique(),
        on="metabolite",
        how="inner",
    )[["sample", "metabolite", "isotopologue", "natural_fraction"]]
    return msts.join(
        natural_with_samples,
        on=["sample", "metabolite", "isotopologue"],
        how="right",
    ).with_columns(
        measured_fraction=close(
            pl.col("measurement"),
            over=["sample", "metabolite"],
        )
    )


def main():
    name_to_func = {"daria": prepare_data_daria, "sergi": prepare_data_sergi}
    for name, raw_files in RAW_FILES.items():
        prepfunc = name_to_func[name]
        raw = (pl.read_csv(rf) for rf in raw_files)
        prepared = prepfunc(*raw)
        prepared.write_csv(OUT_DIR / f"measurements-{name}.csv")


if __name__ == "__main__":
    main()
