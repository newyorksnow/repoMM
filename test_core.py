"""Command-line interface: python -m couplingdx ..."""
import argparse
from .core import design_check, NU_USABLE

def main(argv=None):
    p = argparse.ArgumentParser(
        prog="couplingdx",
        description="Feasibility diagnostic for two-person coupling designs.")
    p.add_argument("--baseline-se", type=float, required=True,
                   help="SE of one individual's baseline from one session")
    p.add_argument("--baseline-contrast", type=float, required=True,
                   help="expected |difference| between the two baselines")
    p.add_argument("--gap", type=float, required=True,
                   help="expected joint-session gap under no coupling")
    p.add_argument("--sigma", type=float, default=None,
                   help="gap measurement error (default: counting floor)")
    p.add_argument("--n-events-total", type=float, default=None,
                   help="total events across both individuals")
    p.add_argument("--dispersion", type=float, default=1.0)
    p.add_argument("--n-baseline-sessions", type=int, default=1)
    p.add_argument("--scale-points", type=int, default=7)
    p.add_argument("--baseline-offset", type=float, default=0.0)
    p.add_argument("--nu-star", type=float, default=NU_USABLE)
    a = p.parse_args(argv)
    print(design_check(
        baseline_se=a.baseline_se, baseline_contrast=a.baseline_contrast,
        gap=a.gap, sigma=a.sigma, n_events_total=a.n_events_total,
        dispersion=a.dispersion, n_baseline_sessions=a.n_baseline_sessions,
        scale_points=a.scale_points, baseline_offset=a.baseline_offset,
        nu_star=a.nu_star))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
