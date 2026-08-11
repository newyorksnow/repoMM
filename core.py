"""Feasibility diagnostic for two-person coupling designs.

A coupling design observes two individuals together and compares the joint
observation against baselines estimated separately for each of them.  The
estimator of the coupling then divides by a difference of two noisy baselines.
This module computes, before data are collected, whether that estimator can
resolve the coupling at all.

Three quantities decide it:

  nu     solo contrast ratio: the separation between the two individuals'
         baselines, in units of the standard error of a baseline estimate.
  gap    the joint-observation gap expected under no coupling, in the units the
         observable is counted in.
  sigma  the measurement error of that gap.  For a counting observable this has
         an irreducible floor at sqrt(total events).

The contraction factor recoverable from a single joint session is

    h(theta) = (1 - exp(-theta)) / theta,      theta = 2 * kappa * T,

and a usable window in theta exists only when  sigma < gap / nu_star.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

__all__ = [
    "h", "inv_h", "ceiling", "counting_floor", "required_duration",
    "design_check", "DesignReport",
]

# Simulation-derived defaults (see reference).  NU_MOMENTS is where the ratio
# estimator's moments begin to exist; NU_USABLE is where recovery is actually
# acceptable.  They differ by a factor of four, and using the former as a design
# criterion is the most common way to get this wrong.
NU_MOMENTS = 3.0
NU_USABLE = 12.0


def h(theta: float) -> float:
    """Contraction factor. h(0)=1, strictly decreasing, ~1/theta for large theta."""
    if theta < 0:
        raise ValueError("theta must be non-negative")
    if theta < 1e-8:
        return 1.0 - theta / 2.0
    return -math.expm1(-theta) / theta


def inv_h(value: float) -> float:
    """Invert h. Returns 0.0 if value >= 1 (no coupling resolvable), inf if <= 0."""
    if value >= 1.0:
        return 0.0
    if value <= 0.0:
        return math.inf
    lo, hi = 1e-12, 1e4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if h(mid) > value:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ceiling(scale_points: int = 7, baseline_offset: float = 0.0) -> float:
    """Largest contrast in h obtainable across an ordinal scale.

    Under a proportional map kappa = c * (baseline_offset + s), the effective
    ratio is q = (offset + max) / (offset + min), and the attainable spread is
    bounded by (sqrt(q) - 1) / (sqrt(q) + 1).  A nonzero offset -- individuals at
    the bottom of the scale still couple somewhat -- lowers the bound sharply.
    """
    if scale_points < 2:
        raise ValueError("scale_points must be at least 2")
    if baseline_offset < 0:
        raise ValueError("baseline_offset must be non-negative")
    q = (baseline_offset + scale_points) / (baseline_offset + 1.0)
    return (math.sqrt(q) - 1.0) / (math.sqrt(q) + 1.0)


def counting_floor(n_events_total: float) -> float:
    """Irreducible error on a difference of two independent counts.

    Holds under perfect measurement and a perfectly specified model; any further
    error source adds to it.  Overdispersion inflates it by sqrt(the dispersion
    ratio), so pass n_events_total * dispersion if events cluster.
    """
    if n_events_total < 0:
        raise ValueError("n_events_total must be non-negative")
    return math.sqrt(n_events_total)


def required_duration(baseline_contrast: float, event_rate_total: float,
                      dispersion: float = 1.0, nu_star: float = NU_USABLE) -> float:
    """Session duration at which a counting design becomes resolvable.

    For a counting observable the expected gap grows as `baseline_contrast * T`
    while its error grows only as sqrt(event_rate_total * dispersion * T), so the
    ratio grows as sqrt(T) and a sufficient duration always exists in principle.
    Setting that ratio to `nu_star` gives

        T_req = nu_star**2 * event_rate_total * dispersion / baseline_contrast**2

    The same expression is the baseline duration per individual needed to reach
    nu_star, because both quantities scale identically in T.  Whether a design is
    practical is therefore decided by contrast relative to the square root of the
    event rate, not by either alone.
    """
    if baseline_contrast <= 0:
        raise ValueError("baseline_contrast must be positive")
    if event_rate_total <= 0:
        raise ValueError("event_rate_total must be positive")
    return nu_star ** 2 * event_rate_total * dispersion / baseline_contrast ** 2


@dataclass
class DesignReport:
    """Result of a feasibility check. `feasible` is the headline."""

    nu: float
    gap: float
    sigma: float
    sigma_is_floor: bool
    nu_star: float
    theta_max: float
    theta_min: float
    max_spread: float
    feasible: bool
    notes: list = field(default_factory=list)

    @property
    def window(self) -> tuple:
        return (self.theta_min, self.theta_max)

    def __str__(self) -> str:
        L = ["Coupling design feasibility", "=" * 34,
             f"  solo contrast ratio  nu   = {self.nu:8.2f}   (need >= {self.nu_star:g})",
             f"  expected gap              = {self.gap:8.2f}",
             f"  gap measurement error     = {self.sigma:8.2f}"
             + ("  [counting floor]" if self.sigma_is_floor else ""),
             f"  max spread across scale   = {self.max_spread:8.3f}", ""]
        if self.feasible:
            L.append(f"  FEASIBLE: usable theta in "
                     f"[{self.theta_min:.3f}, {self.theta_max:.3f}]")
        else:
            L.append("  NOT FEASIBLE: the usable window in theta is empty.")
        if self.notes:
            L += ["", "Notes:"] + [f"  - {n}" for n in self.notes]
        return "\n".join(L)


def design_check(
    baseline_se: float,
    baseline_contrast: float,
    gap: float | None = None,
    sigma: float | None = None,
    n_events_total: float | None = None,
    dispersion: float = 1.0,
    n_baseline_sessions: int = 1,
    scale_points: int = 7,
    baseline_offset: float = 0.0,
    nu_star: float = NU_USABLE,
) -> DesignReport:
    """Check whether a coupling design can resolve the coupling.

    Parameters
    ----------
    baseline_se
        Standard error of one individual's baseline from ONE session.
    baseline_contrast
        Expected |difference| between the two individuals' baselines, in the
        same units.  For random pairing from a population with between-person SD
        s, the median is about 0.95 * s; pre-screening for contrast raises it.
    gap
        Expected joint-observation gap under no coupling.  If omitted it is
        taken as baseline_contrast * duration, which requires `duration`
        semantics the caller supplies by scaling baseline_contrast beforehand.
    sigma
        Gap measurement error.  If omitted, derived from `n_events_total` via
        the counting floor, which is a lower bound -- so an omitted sigma yields
        an optimistic verdict.
    n_events_total
        Total events across both individuals in the joint session.
    dispersion
        Variance-to-mean ratio of the counting process; > 1 if events cluster.
    n_baseline_sessions
        Number of baseline sessions per individual; baseline_se scales as
        1/sqrt(this).
    scale_points, baseline_offset
        Grouping scale, passed to `ceiling`.
    nu_star
        Contrast-ratio threshold.  Defaults to the simulation-derived usable
        value, not the weaker threshold at which moments exist.
    """
    if baseline_se <= 0:
        raise ValueError("baseline_se must be positive")
    if baseline_contrast < 0:
        raise ValueError("baseline_contrast must be non-negative")
    if n_baseline_sessions < 1:
        raise ValueError("n_baseline_sessions must be at least 1")

    notes = []
    se = baseline_se / math.sqrt(n_baseline_sessions)
    sigma_pair = math.sqrt(2.0) * se
    nu = baseline_contrast / sigma_pair

    if gap is None:
        raise ValueError("gap must be supplied")
    if gap <= 0:
        raise ValueError("gap must be positive; a dyad with no expected gap "
                         "carries no information about coupling")

    sigma_is_floor = False
    if sigma is None:
        if n_events_total is None:
            raise ValueError("supply either sigma or n_events_total")
        sigma = counting_floor(n_events_total * dispersion)
        sigma_is_floor = True
        notes.append("sigma taken as the counting floor, a lower bound; the "
                     "verdict is therefore optimistic.")
    if sigma <= 0:
        raise ValueError("sigma must be positive")

    if n_events_total is not None:
        floor = counting_floor(n_events_total * dispersion)
        if sigma < floor - 1e-9:
            notes.append(
                f"supplied sigma ({sigma:.2f}) is below the counting floor "
                f"({floor:.2f}) and is not attainable; using the floor.")
            sigma, sigma_is_floor = floor, True

    theta_max = inv_h(nu_star * sigma / gap)
    # lower end of the noiseless window: half the maximal spread across the scale
    theta_min = 0.115
    feasible = (nu >= nu_star) and (theta_max > theta_min)

    if nu < nu_star:
        need = (nu_star / nu) ** 2 * n_baseline_sessions
        notes.append(f"nu is {nu:.2f}, below {nu_star:g}. Reaching it by "
                     f"replication alone needs about {need:.0f} baseline "
                     f"sessions per individual; raising contrast is usually "
                     f"cheaper and also raises the gap.")
    if theta_max <= theta_min:
        need_gap = nu_star * sigma / h(theta_min)
        notes.append(f"gap of {gap:.1f} is too small against sigma={sigma:.1f}; "
                     f"a non-empty window needs a gap above {need_gap:.1f}.")

    return DesignReport(
        nu=nu, gap=gap, sigma=sigma, sigma_is_floor=sigma_is_floor,
        nu_star=nu_star, theta_max=theta_max, theta_min=theta_min,
        max_spread=ceiling(scale_points, baseline_offset),
        feasible=feasible, notes=notes,
    )
