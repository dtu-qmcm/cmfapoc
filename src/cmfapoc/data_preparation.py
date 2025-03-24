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
    "sergi": RAW_DIR / "PivotTable_freqs 1.csv",
    "daria": RAW_DIR / "processed_calibration_measurements_HEK_LCMS.csv",
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


def prepare_data_daria(raw: pl.DataFrame) -> pl.DataFrame:
    """Prepare Daria's calibration data."""
    outcols = [
        "sample",
        "metabolite",
        "isotopologue",
        "measurement",
        "measured_fraction",
        "natural_fraction",
    ]

    out = raw.rename(
        {
            "met": "metabolite",
            "Theoretical": "theoretical_fraction",
            "Height": "measurement",
            "sample_group": "sample",
        }
    ).with_columns(isotopologue=pl.col("ID").str.tail(1).cast(pl.Int8))
    natural = (
        out.group_by(["metabolite", "isotopologue"])
        .agg(pl.col("theoretical_fraction").first())
        .with_columns(
            natural_fraction=close(pl.col("theoretical_fraction"), "metabolite")
        )
    )
    natural_with_samples = natural.join(
        out[["metabolite", "sample"]].unique(),
        on="metabolite",
        how="inner",
    )
    out = out.join(
        natural_with_samples,
        on=["metabolite", "isotopologue", "sample"],
        how="right",
    )
    out = (
        out.with_columns(
            measured_fraction=close(
                pl.col("measurement"), over=["metabolite", "sample"]
            )
        )
        .select(outcols)
        .sort(["sample", "metabolite", "isotopologue"])
    )
    return out


def main():
    for (name, raw_file), prepfunc in zip(
        RAW_FILES.items(), [prepare_data_sergi, prepare_data_daria]
    ):
        raw = pl.read_csv(raw_file)
        prepared = prepfunc(raw)
        prepared.write_csv(OUT_DIR / f"measurements-{name}.csv")


if __name__ == "__main__":
    main()
