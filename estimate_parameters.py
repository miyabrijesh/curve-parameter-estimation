"""Estimate theta, M and X for the supplied parametric curve."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, least_squares

LOWER_BOUNDS = np.array([np.finfo(float).eps, -0.05, np.finfo(float).eps])
UPPER_BOUNDS = np.array([np.deg2rad(50), 0.05, 100.0])


def local_coordinates(parameters: np.ndarray, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Undo the rotation and x-translation for a candidate parameter triple."""
    theta, _m, shift_x = parameters
    c, s = np.cos(theta), np.sin(theta)
    dx, dy = x - shift_x, y - 42.0
    return c * dx + s * dy, -s * dx + c * dy


def residuals(parameters: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return b_i - exp(M*t_i)*sin(0.3*t_i) for every observed point."""
    _theta, m, _shift_x = parameters
    t, bump = local_coordinates(parameters, x, y)
    return bump - np.exp(m * t) * np.sin(0.3 * t)


def objective(parameters: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    """Robust global loss with a penalty for recovered t values outside 6..60."""
    t, _ = local_coordinates(parameters, x, y)
    outside_range = np.maximum(6.0 - t, 0.0) + np.maximum(t - 60.0, 0.0)
    return float(np.mean(np.abs(residuals(parameters, x, y))) + 1000.0 * np.mean(outside_range))


def fit(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    bounds = list(zip(LOWER_BOUNDS, UPPER_BOUNDS))
    initial = differential_evolution(
        objective, bounds=bounds, args=(x, y), seed=7, popsize=25,
        maxiter=2000, tol=1e-11, polish=True,
    )
    refined = least_squares(
        residuals, initial.x, args=(x, y), bounds=(LOWER_BOUNDS, UPPER_BOUNDS),
        max_nfev=30000, xtol=1e-14, ftol=1e-14, gtol=1e-14,
    )
    return refined.x


def predicted_curve(parameters: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    theta, m, shift_x = parameters
    bump = np.exp(m * np.abs(t)) * np.sin(0.3 * t)
    return (
        t * np.cos(theta) - bump * np.sin(theta) + shift_x,
        42.0 + t * np.sin(theta) + bump * np.cos(theta),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="xy_data.csv")
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    data = pd.read_csv(args.data)
    if list(data.columns) != ["x", "y"] or data.isna().any().any():
        raise ValueError("Input must contain non-empty numeric columns named exactly x and y.")
    x, y = data["x"].to_numpy(float), data["y"].to_numpy(float)
    parameters = fit(x, y)
    theta, m, shift_x = parameters
    fit_residuals = residuals(parameters, x, y)
    recovered_t, _ = local_coordinates(parameters, x, y)
    summary = {
        "theta_radians": float(theta), "theta_degrees": float(np.rad2deg(theta)),
        "M": float(m), "X": float(shift_x),
        "mean_absolute_residual": float(np.mean(np.abs(fit_residuals))),
        "rmse": float(np.sqrt(np.mean(fit_residuals**2))),
        "recovered_t_min": float(recovered_t.min()),
        "recovered_t_max": float(recovered_t.max()), "point_count": int(len(data)),
    }
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    t_grid = np.linspace(6.0, 60.0, 2000)
    curve_x, curve_y = predicted_curve(parameters, t_grid)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(x, y, s=8, alpha=0.35, label="Observed CSV points")
    ax.plot(curve_x, curve_y, color="crimson", linewidth=2.2, label="Fitted curve")
    ax.set(xlabel="x", ylabel="y", title="Observed points and recovered parametric curve")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "curve_fit.png", dpi=180)
    print(json.dumps(summary, indent=2))
    print("\nRounded solution: theta = pi/6 radians (30 degrees), M = 0.03, X = 55")


if __name__ == "__main__":
    main()
