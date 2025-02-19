"""Code for preparing data."""

from collections.abc import Sequence
from pathlib import Path

from polars._typing import IntoExpr
import polars as pl

from cmfapoc.simulation import simulate

ROOT = Path(__file__).parent.parent.parent
RAW_FILE = ROOT / "data" / "raw" / "PivotTable_freqs 1.csv"
OUT_FILE = ROOT / "data" / "prepared" / "measurements.csv"
SIM_ERROR_SD = 0.1


def normalise(expr: pl.Expr, over: IntoExpr | Sequence[IntoExpr]):
    """Normalise expression by sum over given columns.

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
    >>> df.with_columns(c=normalise(pl.col("a"), over=["b"]))
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


def prepare_data(raw: pl.DataFrame) -> pl.DataFrame:
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
            fraction=normalise(
                pl.col("raw_fraction"),
                over=["measurement_id", "species"],
            ),
        )
        .with_columns(is_c12=pl.col("sample_id") == 0)
    )
    return out


def main():
    raw = pl.read_csv(RAW_FILE)
    data = prepare_data(raw)
    for transformation in ("alr", "clr", "ilr"):
        sim = simulate(
            dataset=data,
            transformation=transformation,
            error_sd=SIM_ERROR_SD,
        )
        data = data.with_columns(sim.alias(f"sim_fraction_{transformation}"))
    data.write_csv(OUT_FILE)


if __name__ == "__main__":
    main()
