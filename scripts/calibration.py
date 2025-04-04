import marimo

__generated_with = "0.11.26"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Compositional fluxomics example: calibration

        This notebook presents a simple example of compositional and non compositional analysis of some isotope labelling data used in fluxomics.

        The example relates to calibration, i.e. conducting some measurements of known compositions.

        Specifically, an experiment analysed 18 samples of biological material, resulting in 127 sets of measurements. Each set of measurements contains a mass spectrometer reading for each isotopic mass class, aka "isotopologue", of a metabolite. These sets of measurements are compositional, as the readings carry relative information, but do not indicate the absolute abundances of the measured isotopologues.

        We want to know how accurate our mass spectrometer readings were. Luckily we have a textbook that tells us in what proportion each metabolite's isotopologues occur naturally, so we can compare each set of measurements with what the textbook says they should be.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Introduction

        First we load some libraries and our calibration dataset, then check that the number of measurement sets is as expected. 
        """
    )
    return


@app.cell
def _():
    from cycler import cycler
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import polars as pl

    from matplotlib import pyplot as plt
    from numpy.typing import NDArray
    from skbio.stats.composition import clr, clr_inv, ilr, ilr_inv, alr, alr_inv, sbp_basis

    from cmfapoc.data_preparation import close

    plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('tab10').colors)

    return (
        NDArray,
        Path,
        alr,
        alr_inv,
        close,
        clr,
        clr_inv,
        cycler,
        ilr,
        ilr_inv,
        mo,
        np,
        pl,
        plt,
        sbp_basis,
    )


@app.cell
def _(Path, __file__, pl):
    DATA_FILE = Path(__file__).parent.parent / "data" / "prepared" / "measurements-daria.csv"

    measurements = pl.read_csv(DATA_FILE)
    measurements
    return DATA_FILE, measurements


@app.cell
def _(measurements, pl):
    n_measurement = len(measurements.filter(pl.col("measurement").is_not_null()).group_by(["sample", "metabolite"]).mean())
    print(f"number of measurement sets: {n_measurement}")
    return (n_measurement,)


@app.cell
def _(mo):
    mo.md(
        r"""
        To illustrate what this data tends to look like, the cell below shows an example set of measurements, for metabolite [ribulose 5 phosphate](https://en.wikipedia.org/wiki/Ribulose_5-phosphate) aka "ru5p" in sample "HEK_Wt_QC1_1x_split2_inj1". 

        The measured value for the isotopologue "m0", i.e. unlabelled ru5p, was 134443.22. The measurement for isotopologue "m1" was 7651.59. Closing these measurements over all ru5p isotopologues shows a ratio of approximately 0.946:0.054. This agrees pretty well with the natural fractions, which come from a textbook and suggest an expected ratio of 0.944:0.056.
        """
    )
    return


@app.cell
def _(measurements, pl):
    measurements.filter(
        (pl.col("metabolite") == "ru5p") 
        & (pl.col("sample") == "HEK_Wt_QC1_1x_split2_inj1")
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        There are also some measurement sets where we didn't manage to detect all the isotopologues that should theoretically exist.

        In the example below, isotopologues "m0", "m1", "m2" and "m4" of fructose-1,6-diphosphate were measured, but isotopologue "m3" was not.
        """
    )
    return


@app.cell
def _(measurements, pl):
    measurements.filter(
        (pl.col("metabolite") == "fdp") 
        & (pl.col("sample") == "HEK_Wt_QC1_1x_split1_inj1")
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Quick exploration 

        To get a quick idea about the agreement between  we can start by just dropping all the null rows and re-closing the measurement columns.
        """
    )
    return


@app.cell
def _(close, measurements, pl):
    msts_no_zeros = (
        measurements
        .filter(
            (pl.col("measurement").is_not_null())
            & (pl.col("metabolite") != "akg")
            & (pl.col("measurement").count().over("sample", "metabolite") > 1))
        .with_columns(
            natural_fraction=close(pl.col("natural_fraction"), over=("sample", "metabolite")),
            measured_fraction=close(pl.col("measurement"), over=("sample", "metabolite")),
        )
    )
    return (msts_no_zeros,)


@app.cell
def _(mo):
    mo.md(
        r"""
        The next step is to transform the measurements by centered log ratio, additive log ratio and isometric log ratio, and then look for patterns in the differences between observed and expected measurements.

        We are particularly interested in patterns in the errors that depend on the absolute value of the measurement. Mass spectrometers can produce unexpected measurements for very high or very low target values. At the high end, the measurement is likely to be truncated above at the spectrometer's maximum detectable value. At the low end, the noise can increase as the target value approaches the lower detection limit. There can also be measurement-value-dependent biases!
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ### Centered log ratio

        Applying the CLR transformation to the non-zero measurements shows an interesting pattern: the residuals are more dispersed when the absolute measurement is lower. This is roughly what we expect given the spectrometer's lower detection limit.

        However, it is hard to diagnose biases from the CLR data because it is centered around zero per measurement, as shown in the cell below. 
        """
    )
    return


@app.cell
def _(clr, msts_no_zeros, pl):
    def func_clr(subdf: pl.DataFrame) -> pl.DataFrame:
        subdf = subdf.sort("isotopologue")
        row0 = subdf.row(0, named=True)
        met = row0["metabolite"]
        sample = row0["sample"]
        total_mst = subdf["measurement"]
        ratio = [f"m{i}:gmean" for i in range(len(subdf))]
        m = clr(subdf["measurement"])
        nf = clr(subdf["natural_fraction"])
        resid = m - nf
        return pl.DataFrame(
            {
                "metabolite": met,
                "sample": sample,
                "ratio": ratio,
                "m": m,
                "nf": nf,
                "resid": resid,
                "total_measurement": total_mst,
            }
        )



    msts_no_zeros_clr = (
        msts_no_zeros
        .group_by("sample", "metabolite", maintain_order=True)
        .map_groups(func_clr)
    )
    msts_no_zeros_clr.filter(metabolite="fdp", sample="HEK_Wt_QC1_1x_split1_inj1")
    return func_clr, msts_no_zeros_clr


@app.cell
def _(msts_no_zeros_clr, plt):
    def plot_sct(df, title):
        f, ax = plt.subplots()
        for ((met,), subdf) in df.group_by("ratio", maintain_order=True):
            ax.scatter(subdf["total_measurement"], subdf["resid"], s=10, label=met)
        f.legend(loc="right", bbox_to_anchor=(1.15, 0.5), title="Ratio")
        ax.semilogx()
        ax.set(xlabel="Total measurement (raw units)", ylabel="Residual", title=title)

    plot_sct(msts_no_zeros_clr, "CLR residuals")
    plt.show()
    return (plot_sct,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ### Additive Log Ratio

        Next we do the same thing using the additive log ratio transformation, i.e. finding the log ratio of each labelled isotopologue with  the unlabelled isotopologue "m0".

        This results in a pretty inconclusive plot - we're not really sure why! Maybe it's because the ALR transformation is not isometric?
        """
    )
    return


@app.cell
def _(alr, msts_no_zeros, pl):
    def func_alr(subdf: pl.DataFrame) -> pl.DataFrame:
        subdf = subdf.sort("isotopologue")
        row0 = subdf.row(0, named=True)
        met = row0["metabolite"]
        sample = row0["sample"]
        total_mst = row0["measurement"] + subdf["measurement"].tail(len(subdf)-1)
        ratio = [f"m{i+1}:m0" for i in range(len(subdf)-1)]
        m = alr(subdf["measurement"], denominator_idx=0)
        nf = alr(subdf["natural_fraction"], denominator_idx=0)
        resid = m - nf
        return pl.DataFrame(
            {
                "metabolite": met,
                "sample": sample,
                "ratio": ratio,
                "m": m,
                "nf": nf,
                "resid": resid,
                "total_measurement": total_mst,
            }
        )

    msts_no_zeros_alr = (
        msts_no_zeros
        .group_by("sample", "metabolite", maintain_order=True)
        .map_groups(func_alr)
    )
    msts_no_zeros_alr
    return func_alr, msts_no_zeros_alr


@app.cell
def _(msts_no_zeros_alr, plot_sct, plt):
    plot_sct(msts_no_zeros_alr, "ALR residuals")
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ### Isometric Log Ratio

        Finally, here is the same plot for the data transformed using isometric log ratios.

        To define an orthonormal basis we use sequential binary partitioning, starting with "m0" vs all labelled isotopologues, then "m1" vs all isotopologues with more than one label, and so on.
        """
    )
    return


@app.cell
def _(NDArray, ilr, msts_no_zeros, np, pl, sbp_basis):
    def get_orthonormal_basis_and_labels(size: int) -> tuple[NDArray, list[str]]:
        if size < 2:
            raise ValueError("size must be greater than 1")
        out = np.zeros((size-1, size))
        labels = [f"m{i+1}+:m{i}" for i in range(size-1)]
        for row_ix in range(size-1):
            out[row_ix, row_ix] = -1
            out[row_ix, row_ix+1:] = 1
        return sbp_basis(out), labels


    def func_ilr(subdf: pl.DataFrame) -> pl.DataFrame:
        subdf = subdf.sort("isotopologue")
        row0 = subdf.row(0, named=True)
        met = row0["metabolite"]
        sample = row0["sample"]
        total_mst = subdf["measurement"].cum_sum(reverse=True).head(len(subdf)-1)
        basis, ratios = get_orthonormal_basis_and_labels(len(subdf))
        m = ilr(subdf["measurement"], basis=basis)
        nf = ilr(subdf["natural_fraction"], basis=basis)
        resid = m - nf
        return pl.DataFrame(
            {
                "metabolite": met,
                "sample": sample,
                "ratio": ratios,
                "m": m,
                "nf": nf,
                "resid": resid,
                "total_measurement": total_mst,
            }
        )

    msts_no_zeros_ilr = (
        msts_no_zeros
        .group_by("sample", "metabolite", maintain_order=True)
        .map_groups(func_ilr)
    )
    msts_no_zeros_ilr
    return func_ilr, get_orthonormal_basis_and_labels, msts_no_zeros_ilr


@app.cell
def _(msts_no_zeros_ilr, plot_sct, plt):
    plot_sct(msts_no_zeros_ilr, "ILR residuals")
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""This produces a very interesting graph! Ratios with total absolute measurements more than about $10^5$ seem to be mor or less unbiased, whereas lower-valued ratios tend to be systematically biased towards the more labelled side of the ratio.""")
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## What next?

        We aren't totally sure how we should model the apparent bias in the ILR measurements. One option is to include the total measurement as a covariate in a model of the transformed measurements. Another would be to make a model with both absolute and compositional modules. The second option would potentially allow a more explicit model of the biasing process, but seems like it might be a bit trickier to implement. 
        """
    )
    return


if __name__ == "__main__":
    app.run()
