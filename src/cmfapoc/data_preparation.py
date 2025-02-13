"""Code for preparing data."""

from pathlib import Path
import polars as pl

ROOT = Path(__file__).parent.parent.parent
RAW_FILE = ROOT / "data" / "raw" / "PivotTable_freqs 1.csv"
OUT_FILE = ROOT / "data" / "prepared" / "measurements.csv"


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
            fraction=pl.col("raw_fraction")
            / (
                pl.col("raw_fraction").sum().over(["measurement_id", "species"])
            ),
        )
        .with_columns(is_c12=pl.col("sample_id") == 0)
    )
    return out


if __name__ == "__main__":
    raw = pl.read_csv(RAW_FILE)
    data = prepare_data(raw)
    data.write_csv(OUT_FILE)
