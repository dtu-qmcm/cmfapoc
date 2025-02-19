"""Code for simulating fake compositional data, given a real dataset.

The simulation is done by transforming the measurements to unconstrained space using either the centered log ratio (clr), the additive log ratio (alr), or the isometric log ratio (ilr) transformations, then adding Gaussian noise to the transformed data, and finally transforming the data back to compositional space using the inverse of the transformation.

The table containing the dataset must have the following columns:
- fraction: the component of the compositional measurement for this row.
- measurement_id: the id of the measurement.
- species: the name of the species.

"""

from functools import partial
import numpy as np
import polars as pl
from skbio.stats.composition import alr, alr_inv, clr, clr_inv, ilr, ilr_inv

TRANSFORMATIONS = {
    "alr": (alr, alr_inv),
    "clr": (clr, clr_inv),
    "ilr": (ilr, ilr_inv),
}


def simulate_compositional_measurement(
    df: pl.DataFrame,
    error_sd: float,
    transformation: str = "clr",
) -> pl.DataFrame:
    """Simulate a single compositional measurement.

    Args:
        df: A dataframe with a single measurement in the column "fraction".
        error_sd: The standard deviation of the error, on centered-log-ratio scale.
        transformation: name of the transformation to use. Options are "alr", "clr", and "ilr".

    Returns:
        A series with the simulated measurements.
    """
    if transformation not in TRANSFORMATIONS:
        options = str(list(TRANSFORMATIONS.keys()))
        msg = f"Parameter 'transformation' must be one of: {options}."
        raise ValueError(msg)
    f, finv = TRANSFORMATIONS[transformation]
    unconstrained = f(df["fraction"])
    deviations = np.random.normal(size=unconstrained.shape) * error_sd
    sim_unconstrained = unconstrained + deviations
    return pl.DataFrame({"sim_fraction": finv(sim_unconstrained)})


def simulate(
    dataset: pl.DataFrame,
    error_sd: float = 0.1,
    transformation: str = "clr",
) -> pl.Series:
    """Apply simulate_compositional_measurement to a dataset.

    Args:
        dataset: A dataframe with columns "fraction", "measurement_id" and "species".
        error_sd: The standard deviation of the error, on centered-log-ratio scale.
        transformation: name of the transformation to use. Options are "alr", "clr", and "ilr".

    Returns:
        A series with the simulated measurements.
    """
    sim_func = partial(
        simulate_compositional_measurement,
        error_sd=error_sd,
        transformation=transformation,
    )
    return dataset.group_by(
        "measurement_id",
        "species",
        maintain_order=True,
    ).map_groups(sim_func)["sim_fraction"]


if __name__ == "__main__":
    from pathlib import Path

    root = Path(__file__).parent.parent.parent
    df = pl.read_csv(root / "data" / "prepared" / "measurements.csv")
    sim = simulate(df, transformation="alr", error_sd=0.5)
    print(df.with_columns(simulated_measurement=sim))
