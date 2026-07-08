#!/usr/bin/env python3
"""
Analyze A/B test results from a CSV of events.

Expected columns:
  variant     - "A" or "B" (any two labels work; the first seen becomes control)
  converted   - 0 or 1

Usage:
  python3 analyze.py sample_data.csv
  python3 analyze.py sample_data.csv --alpha 0.01
  python3 analyze.py sample_data.csv --bayes-samples 50000

No third-party dependencies. Runs on stdlib Python 3.8+.
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from collections import defaultdict


# ---- Normal distribution (same primitives as stats.js) ---------------------

def norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def norm_inv(p: float) -> float:
    """Beasley-Springer-Moro inverse normal CDF."""
    if not 0.0 < p < 1.0:
        raise ValueError("norm_inv: p must be in (0, 1)")
    a = [-39.6968302866538, 220.946098424521, -275.928510446969,
         138.357751867269, -30.6647980661472, 2.50662827745924]
    b = [-54.4760987982241, 161.585836858041, -155.698979859887,
         66.8013118877197, -13.2806815528857]
    c = [-0.00778489400243029, -0.322396458041136, -2.40075827716184,
         -2.54973253934373, 4.37466414146497, 2.93816398269878]
    d = [0.00778469570904146, 0.32246712907004, 2.445134137143,
         3.75440866190742]
    plow = 0.02425
    phigh = 1 - plow
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return ((((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1))
    if p <= phigh:
        q = p - 0.5
        r = q * q
        return (((((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q) /
                (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1))
    q = math.sqrt(-2 * math.log(1 - p))
    return -((((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1))


# ---- Analyses --------------------------------------------------------------

def chi_sq_survival_1(x: float) -> float:
    """P(χ²_1 > x) via the identity χ²_1 = Z²."""
    if x <= 0:
        return 1.0
    return 2 * (1 - norm_cdf(math.sqrt(x)))


def srm_check(nA: int, nB: int, expected_a: float = 0.5) -> dict:
    """Sample Ratio Mismatch chi-square test. Alerts at p < 0.001 (Kohavi)."""
    total = nA + nB
    if total == 0:
        return {"chi_sq": 0.0, "p_value": 1.0, "alert": False}
    eA = total * expected_a
    eB = total * (1 - expected_a)
    chi = (nA - eA) ** 2 / eA + (nB - eB) ** 2 / eB
    p = chi_sq_survival_1(chi)
    return {"chi_sq": chi, "p_value": p, "alert": p < 0.001,
            "observed_a": nA / total, "expected_a": expected_a}


def always_valid(nA: int, xA: int, nB: int, xB: int,
                 tau: float = 0.02, alpha: float = 0.05) -> dict:
    """mSPRT with Gaussian mixing prior N(0, τ²) on the true lift.
    Johari, Koomen, Pekelis, Walsh (2015) — always-valid inference."""
    pA = xA / nA
    pB = xB / nB
    delta = pB - pA
    V = pA * (1 - pA) / nA + pB * (1 - pB) / nB
    if V == 0:
        return {"p_value": 1.0, "ci": (delta, delta), "half_width": 0.0, "tau": tau}
    t2 = tau * tau
    log_lambda = 0.5 * math.log(V / (V + t2)) + (delta * delta * t2) / (2 * V * (V + t2))
    p_value = min(1.0, math.exp(-log_lambda))
    half_width = math.sqrt(
        (2 * V * (V + t2) / t2) * (math.log(1 / alpha) + 0.5 * math.log((V + t2) / V))
    )
    return {"p_value": p_value, "ci": (delta - half_width, delta + half_width),
            "half_width": half_width, "tau": tau}


def frequentist(nA: int, xA: int, nB: int, xB: int, alpha: float = 0.05) -> dict:
    pA = xA / nA
    pB = xB / nB
    p_pool = (xA + xB) / (nA + nB)
    se_pool = math.sqrt(p_pool * (1 - p_pool) * (1 / nA + 1 / nB))
    z = 0.0 if se_pool == 0 else (pB - pA) / se_pool
    p_value = 2 * (1 - norm_cdf(abs(z)))
    se_diff = math.sqrt(pA * (1 - pA) / nA + pB * (1 - pB) / nB)
    z_alpha = norm_inv(1 - alpha / 2)
    ci = (pB - pA - z_alpha * se_diff, pB - pA + z_alpha * se_diff)
    return {
        "rateA": pA, "rateB": pB,
        "abs_lift": pB - pA,
        "rel_lift": (pB - pA) / pA if pA else float("nan"),
        "z": z, "p_value": p_value, "ci95": ci,
        "significant": p_value < alpha,
    }


def bayesian(nA: int, xA: int, nB: int, xB: int, samples: int = 20000) -> dict:
    aA, bA = 1 + xA, 1 + nA - xA
    aB, bB = 1 + xB, 1 + nB - xB
    wins = 0
    loss_a = 0.0
    loss_b = 0.0
    draws_a = []
    draws_b = []
    for _ in range(samples):
        a = random.betavariate(aA, bA)
        b = random.betavariate(aB, bB)
        draws_a.append(a)
        draws_b.append(b)
        if b > a:
            wins += 1
            loss_a += b - a
        else:
            loss_b += a - b
    draws_a.sort()
    draws_b.sort()

    def q(arr, qq):
        return arr[int(qq * (len(arr) - 1))]

    return {
        "prob_B_beats_A": wins / samples,
        "expected_loss_pick_A": loss_a / samples,
        "expected_loss_pick_B": loss_b / samples,
        "posterior_A_95": (q(draws_a, 0.025), q(draws_a, 0.975)),
        "posterior_B_95": (q(draws_b, 0.025), q(draws_b, 0.975)),
    }


# ---- CSV loading -----------------------------------------------------------

def load_csv(path: str) -> dict:
    totals = defaultdict(lambda: [0, 0])  # variant -> [visitors, conversions]
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if "variant" not in reader.fieldnames or "converted" not in reader.fieldnames:
            sys.exit("CSV must have columns: variant, converted")
        for row in reader:
            v = row["variant"].strip()
            c = int(row["converted"])
            totals[v][0] += 1
            totals[v][1] += c
    if len(totals) != 2:
        sys.exit(f"Expected exactly 2 variants, got {list(totals.keys())}")
    return dict(totals)


# ---- CLI -------------------------------------------------------------------

def pct(x: float, d: int = 3) -> str:
    return f"{x*100:.{d}f}%"


def main() -> None:
    ap = argparse.ArgumentParser(description="Analyze A/B test CSV.")
    ap.add_argument("csv", help="Path to CSV with columns variant,converted")
    ap.add_argument("--alpha", type=float, default=0.05, help="Significance level (default 0.05)")
    ap.add_argument("--expected-a", type=float, default=0.5, help="Expected traffic share for variant A (default 0.5)")
    ap.add_argument("--tau", type=float, default=0.02, help="mSPRT prior scale on absolute lift (default 0.02)")
    ap.add_argument("--bayes-samples", type=int, default=20000, help="Monte Carlo draws for Bayes (default 20000)")
    ap.add_argument("--seed", type=int, default=None, help="RNG seed for reproducible Bayes draws")
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    totals = load_csv(args.csv)
    variants = sorted(totals.keys())  # A, B alphabetically; predictable
    a_label, b_label = variants[0], variants[1]
    nA, xA = totals[a_label]
    nB, xB = totals[b_label]

    print(f"\n== A/B test: {args.csv} ==")
    print(f"  {a_label}: {nA:,} visitors, {xA:,} conversions ({pct(xA/nA)})")
    print(f"  {b_label}: {nB:,} visitors, {xB:,} conversions ({pct(xB/nB)})")

    srm = srm_check(nA, nB, expected_a=args.expected_a)
    print("\n-- Sample Ratio Mismatch (χ²) --")
    print(f"  Observed {a_label} share: {pct(srm['observed_a'], 2)}   expected: {pct(args.expected_a, 1)}")
    print(f"  χ² = {srm['chi_sq']:.3f}   p = {srm['p_value']:.4g}")
    if srm["alert"]:
        print(f"  ALERT: p < 0.001 — likely bug in assignment/logging. Metric results below are untrustworthy.")
    else:
        print(f"  OK: assignment ratio is consistent with the target split.")

    f = frequentist(nA, xA, nB, xB, alpha=args.alpha)
    print("\n-- Frequentist z-test --")
    print(f"  Absolute lift ({b_label} - {a_label}): {pct(f['abs_lift'])}")
    print(f"  Relative lift:                          {pct(f['rel_lift'], 2)}")
    print(f"  z-score:                                {f['z']:+.3f}")
    print(f"  p-value (two-sided):                    {f['p_value']:.4g}")
    print(f"  95% CI on lift:                         [{pct(f['ci95'][0])}, {pct(f['ci95'][1])}]")
    verdict = "SIGNIFICANT" if f["significant"] else "not significant"
    print(f"  Verdict at α={args.alpha}:                    {verdict}")

    av = always_valid(nA, xA, nB, xB, tau=args.tau, alpha=args.alpha)
    print("\n-- Always-valid (mSPRT, Johari et al.) --")
    print(f"  Always-valid p-value:                   {av['p_value']:.4g}")
    print(f"  Confidence sequence on lift:            [{pct(av['ci'][0])}, {pct(av['ci'][1])}]")
    print(f"  Prior scale τ:                          {pct(args.tau, 2)} absolute")
    av_verdict = "SIGNIFICANT" if av["p_value"] < args.alpha else "not significant"
    print(f"  Verdict at α={args.alpha} (peek-safe):        {av_verdict}")

    b = bayesian(nA, xA, nB, xB, samples=args.bayes_samples)
    print("\n-- Bayesian (Beta(1,1) prior) --")
    print(f"  P({b_label} beats {a_label}):                          {pct(b['prob_B_beats_A'], 1)}")
    print(f"  Expected loss picking {a_label}:               {pct(b['expected_loss_pick_A'])}")
    print(f"  Expected loss picking {b_label}:               {pct(b['expected_loss_pick_B'])}")
    print(f"  95% CrI on {a_label} rate:                     [{pct(b['posterior_A_95'][0])}, {pct(b['posterior_A_95'][1])}]")
    print(f"  95% CrI on {b_label} rate:                     [{pct(b['posterior_B_95'][0])}, {pct(b['posterior_B_95'][1])}]")
    print()


if __name__ == "__main__":
    main()
