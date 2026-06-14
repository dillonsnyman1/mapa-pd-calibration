"""Shared plotting helper for the demo scripts."""

import matplotlib.pyplot as plt


def plot_calibration(bands, smoothed_points, title, out_path):
    """Plot an unsmoothed (step) and smoothed score-to-PD curve.

    Args:
        bands: iterable of (score_min, score_max, pd) tuples.
        smoothed_points: iterable of (score, pd) tuples, densely sampled.
        title: plot title.
        out_path: path to save the PNG to.
    """
    step_x, step_y = [], []
    for score_min, score_max, pd in bands:
        step_x += [score_min, score_max]
        step_y += [pd, pd]

    smooth_x = [score for score, _ in smoothed_points]
    smooth_y = [pd for _, pd in smoothed_points]

    plt.figure(figsize=(9, 6))
    plt.step(step_x, step_y, where="post", label="Unsmoothed (pooled bands)", color="tab:blue")
    plt.plot(smooth_x, smooth_y, label="Smoothed (log-odds interpolation)", color="tab:orange")
    plt.xlabel("Score")
    plt.ylabel("PD")
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved plot to {out_path}")
