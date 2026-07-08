# Splitcheck

A small, self-contained toolkit for planning and analyzing A/B tests
correctly. Zero dependencies — the browser app is plain HTML + vanilla ES
modules, the CLI is stdlib Python 3.

Splitcheck is opinionated about two things most tutorials skip:

1. **Check the split before you trust the metric.** Every result view leads
   with a Sample-Ratio-Mismatch chi-square test. If the assignment is broken,
   nothing downstream is trustworthy — Kohavi's rule.
2. **Fixed-horizon and always-valid are different animals.** You get both:
   the ordinary z-test that requires you to pre-commit a sample size, and the
   mSPRT-based always-valid p-value and confidence sequence, which let you
   peek as often as you like without inflating the false-positive rate.

Everything renders in a browser with no build step. The Python CLI reads a
`variant,converted` CSV and prints the same numbers.

---

## Contents

- [Quick start](#quick-start)
- [What each tab does](#what-each-tab-does)
- [The statistics — full derivations](#the-statistics--full-derivations)
  - [1. Two-proportion z-test](#1-two-proportion-z-test)
  - [2. Sample-size formula](#2-sample-size-formula)
  - [3. Bayesian Beta-Binomial with expected loss](#3-bayesian-beta-binomial-with-expected-loss)
  - [4. Sampling Beta variates: Marsaglia-Tsang](#4-sampling-beta-variates-marsaglia-tsang)
  - [5. Always-valid inference: mSPRT](#5-always-valid-inference-msprt)
  - [6. Sample Ratio Mismatch: chi-square](#6-sample-ratio-mismatch-chi-square)
- [Worked example: `sample_data.csv`](#worked-example-sample_datacsv)
- [Design decisions](#design-decisions)
- [What Splitcheck is NOT](#what-splitcheck-is-not)
- [Files](#files)
- [References](#references)

---

## Quick start

Open `index.html` in a browser (double-click) or serve the folder:

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

Run the CLI on the shipped synthetic data (10k / 10k, true lift ~12% rel):

```bash
python3 analyze.py sample_data.csv --seed 1
```

CLI flags: `--alpha`, `--expected-a`, `--tau`, `--bayes-samples`, `--seed`.

## What each tab does

### 1. Plan — required sample size

Given a baseline conversion rate `p₁` and the smallest relative lift you care
about (MDE), compute how many visitors per arm you need. Uses the two-sided
two-proportion power formula (§2 below). Rule of thumb: **halving the MDE
roughly quadruples the required sample.**

### 2. Analyze — verdict on a finished experiment

Enter totals for A and B. You get four analyses stacked:

- **SRM check** (top): chi-square goodness-of-fit against the target split.
  At `p < 0.001` a red banner tells you to stop trusting the metric results
  until you find the assignment bug.
- **Frequentist**: pooled two-proportion z-test for the p-value; unpooled SE
  for the 95% CI on the lift. Only valid if the sample size was pre-committed.
- **Always-valid (mSPRT)**: an anytime-valid p-value and confidence sequence
  from Johari–Koomen–Pekelis–Walsh (2015). Peek as often as you want. Costs a
  wider interval up front.
- **Bayesian**: Beta(1,1) → Beta posteriors, 20,000 Monte Carlo draws for
  `P(B > A)`, expected loss under each choice, 95% credible intervals.

### 3. Peeking simulator

Sets A and B to the *same* true rate, so any declared winner is a false
positive. Simulates three experimenters:

- **Honest**: one test at `n_max`. Hits nominal α.
- **Naive peeker**: fixed-horizon test at every check, stops on first
  `p < α`. At α = 0.05 with 50 looks, ~30% false positives.
- **Always-valid peeker**: mSPRT test at every check. Stays at or below α no
  matter how often it peeks. This is the fix.

---

## The statistics — full derivations

### 1. Two-proportion z-test

**Setup.** Two arms, `n_A` visitors with `x_A` conversions, `n_B` visitors
with `x_B`. Sample proportions:

```
p̂_A = x_A / n_A
p̂_B = x_B / n_B
Δ̂  = p̂_B - p̂_A       (the observed lift)
```

Splitcheck uses **different standard errors for the p-value and the CI**.
This is subtle but standard, and many tutorials mix it up.

**Pooled SE — for the p-value.**
The p-value tests the null `H₀ : p_A = p_B`. *Under H₀* the two proportions
are equal, so the best variance estimate pools them:

```
p̂    = (x_A + x_B) / (n_A + n_B)      (pooled proportion)
SE₀  = √[ p̂(1 - p̂) · (1/n_A + 1/n_B) ]
z    = Δ̂ / SE₀
p    = 2 · (1 - Φ(|z|))                 (two-sided)
```

**Unpooled SE — for the confidence interval.**
The CI does *not* assume the null. It asks "what values of the true lift are
consistent with what we saw?", so we use each arm's own variance:

```
SE_Δ = √[ p̂_A(1 - p̂_A)/n_A  +  p̂_B(1 - p̂_B)/n_B ]
CI   = Δ̂  ±  z_{α/2} · SE_Δ
```

**Why the split.** If you use pooled SE for both, you get a CI that's too
narrow under alternatives; if you use unpooled SE for both, your test can
disagree with your CI on the boundary. Splitting matches how R's
`prop.test`, most stats textbooks, and every serious A/B testing tool
report the numbers.

The implementation: [`stats.js:analyzeFrequentist`](stats.js) /
[`analyze.py:frequentist`](analyze.py).

### 2. Sample-size formula

Two-proportion z-test, two-sided, per-arm sample size `n` such that a true
effect of size `δ = p₂ − p₁` is detected with power `1 − β` at level `α`:

```
                  (z_{α/2} + z_β)² · [ p₁(1 - p₁) + p₂(1 - p₂) ]
    n_per_arm  =  ────────────────────────────────────────────────
                                   (p₂ - p₁)²
```

**Where the terms come from.** Under `H₀ : p₁ = p₂` the standardized test
statistic `z = (p̂₂ - p̂₁) / SE₀` is approximately N(0, 1). Under `H₁` with
true lift `δ`, `z ≈ N(δ/SE, 1)` where `SE² = p₁(1 - p₁)/n + p₂(1 - p₂)/n`.
For a two-sided test at level α, we reject when `|z| > z_{α/2}`. Ignoring
the tiny lower-tail rejection probability (standard approximation when
`δ > 0`), power is

```
1 - β = Φ(δ/SE - z_{α/2})
      ⇒ δ/SE - z_{α/2} = z_β
      ⇒ δ² = (z_{α/2} + z_β)² · SE²
      ⇒ δ² = (z_{α/2} + z_β)² · [ p₁(1-p₁) + p₂(1-p₂) ] / n
```

which rearranges to the boxed formula. Splitcheck uses unpooled variance
here (the standard convention for sample-size planning under the
alternative, not the null).

**Consequence.** `n ∝ 1/δ²`. Halving the MDE quadruples the required sample.
At `p₁ = 5%, MDE = 10% relative (so δ = 0.5pp), α = 0.05, power = 0.8`,
you need **31,231 per arm**. This matches Evan Miller's calculator to the
last digit.

Implementation: [`stats.js:sampleSizePerArm`](stats.js).

### 3. Bayesian Beta-Binomial with expected loss

**Model.** Independent Beta priors on each arm's true conversion rate:

```
p_A ~ Beta(α₀, β₀)     p_B ~ Beta(α₀, β₀)
x_A | p_A ~ Binomial(n_A, p_A)
x_B | p_B ~ Binomial(n_B, p_B)
```

Splitcheck uses `α₀ = β₀ = 1` (flat prior on [0, 1]). Conjugacy gives the
posteriors immediately:

```
p_A | data ~ Beta(1 + x_A, 1 + n_A - x_A)
p_B | data ~ Beta(1 + x_B, 1 + n_B - x_B)
```

**P(B > A).** No closed form for arbitrary parameters, so Splitcheck draws
20,000 Monte Carlo samples from each posterior and counts:

```
P(B > A) ≈ (1/N) · Σᵢ 𝟙{ p_Bᵢ > p_Aᵢ }
```

**Expected loss.** The Stucchio (2015) decision rule at VWO. If the true
lift is `Δ = p_B - p_A`, the loss from picking A is `max(0, Δ)` — you gave
up the good arm — and the loss from picking B is `max(0, -Δ)`. Under the
posterior:

```
E[loss | pick A] = E[ max(0, p_B - p_A) ]
E[loss | pick B] = E[ max(0, p_A - p_B) ]
```

Both are estimated by Monte Carlo. You **stop the experiment when the
minimum of the two expected losses drops below your "threshold of caring"**
— e.g. 0.1 percentage points of conversion. Both losses being small is the
Bayesian analog of "no meaningful difference detected."

Implementation: [`stats.js:analyzeBayesian`](stats.js).

**Honest caveat.** A Bayesian stopping rule based on `P(B > A) > threshold`
is **not** automatically peek-safe (Robinson 2015, Georgiev 2017). The
expected-loss stopping rule with a fixed threshold is better-behaved but
still not formally always-valid without an additional decision-theoretic
argument. For a formally always-valid answer, use the mSPRT (§5).

### 4. Sampling Beta variates: Marsaglia-Tsang

Both frames rely on Beta samples. Splitcheck's JS implementation uses the
Marsaglia–Tsang (2000) rejection sampler for Gamma, then forms Beta as a
ratio.

```
X ~ Gamma(α, 1),  Y ~ Gamma(β, 1)   ⇒   X / (X + Y) ~ Beta(α, β)
```

**Marsaglia-Tsang for `Gamma(α ≥ 1, 1)`.** Let `d = α - 1/3`,
`c = 1/√(9d)`. Repeat:

1. Draw `Z ~ N(0, 1)`, set `v = (1 + cZ)³`.
2. If `v ≤ 0`, restart.
3. Draw `U ~ U(0, 1)`.
4. Accept `dv` if `U < 1 - 0.0331·Z⁴` (the "squeeze"), or if
   `log U < ½Z² + d(1 - v + log v)` (the "full" test).

Very high acceptance rate, no transcendentals in the common path, no
special-casing near α = 1. For `α < 1`, use the Ahrens–Dieter boost:
`Gamma(α) ~ Gamma(α + 1) · U^(1/α)`.

The `N(0, 1)` draw uses Box–Muller from two `U(0, 1)` samples.

The Python CLI just calls `random.betavariate` — same distribution, simpler
code.

Implementation: [`stats.js:sampleGamma`, `sampleBeta`](stats.js).

### 5. Always-valid inference: mSPRT

The most important piece and the one worth reading carefully.

**Motivation.** A fixed-horizon p-value is only valid at the sample size
you pre-committed to. If you check the dashboard mid-experiment and stop
early on a favorable result, your false-positive rate is *not* α. The
peeking simulator shows the empirical inflation (from 5% to ~30% with 50
looks). The fix is a test statistic whose distribution under H₀ is
controlled at *every* stopping time.

**Ville's inequality.** If `(M_n)` is a non-negative supermartingale under
H₀ starting at `M_0 = 1`, then for any `α ∈ (0, 1)`:

```
P₀( sup_n  M_n  ≥  1/α )   ≤   α
```

So if we build a statistic `Λ_n` that is a non-negative martingale under
H₀ and start it at 1, "reject when `Λ_n ≥ 1/α` for any `n`" is a valid
level-α test. The "for any `n`" is what makes it peek-tolerant.

**The mSPRT construction (Johari–Koomen–Pekelis–Walsh, 2015).** We build
`Λ_n` as a mixture likelihood ratio. Choose a prior `π` over the alternative
parameter (the true lift `δ`). Then:

```
Λ_n  =  ∫  [ p_δ(data₁, …, data_n) / p₀(data₁, …, data_n) ]  π(dδ)
```

is a martingale under H₀, because each single-δ likelihood ratio is a
martingale and mixtures of martingales are martingales.

**Gaussian mixing prior.** Splitcheck uses `π = N(0, τ²)`. Treat the
estimator `Δ̂_n` of the lift as approximately normal with variance
`V_n = p̂_A(1 - p̂_A)/n_A + p̂_B(1 - p̂_B)/n_B` (unpooled — this is the
variance under the alternative, where the mixture lives). Then

```
Λ_n(0)  =  ∫  N(Δ̂_n ; δ, V_n) · N(δ ; 0, τ²)  dδ  /  N(Δ̂_n ; 0, V_n)
```

The integral in the numerator is a Gaussian convolution with variance
`V_n + τ²`, so `Δ̂_n ~ N(0, V_n + τ²)` under the mixture. Both the numerator
and denominator are Gaussian densities at `Δ̂_n`, and after algebra:

```
              ┌─────────────────┐              ┌───────────────────────┐
              │      V_n         │              │  Δ̂_n²  ·  τ²           │
    Λ_n(0) = √│ ─────────────── │  · exp       │ ──────────────────── │
              │    V_n + τ²      │              │  2 V_n · (V_n + τ²)   │
              └─────────────────┘              └───────────────────────┘
```

The **always-valid p-value** at time `n` is

```
p*_n  =  min( 1,  1 / Λ_n(0) )
```

Reject `H₀` the first time this falls below α. By Ville's inequality, the
false-positive rate is at most α no matter how often you peek.

**Confidence sequence.** For each candidate `δ`, run the same construction
with the mixing prior centered at `δ`:

```
                ┌────────────────┐           ┌──────────────────────────┐
                │     V_n         │           │  (Δ̂_n - δ)² ·  τ²          │
    Λ_n(δ) = √  │ ──────────── │  · exp     │ ────────────────────── │
                │   V_n + τ²      │           │   2 V_n · (V_n + τ²)    │
                └────────────────┘           └──────────────────────────┘
```

The `1 - α` confidence sequence is `{ δ : Λ_n(δ) < 1/α }`. Solving:

```
                     ┌─────────────────────────────────────────────────────┐
                     │  2 V_n (V_n + τ²)      ┌                    ┐        │
    δ  ∈  Δ̂_n  ± √   │ ───────────────── ·  │  log(1/α) + ½·log(v)│         │
                     │        τ²              └                    ┘        │
                     └─────────────────────────────────────────────────────┘

    where v = (V_n + τ²) / V_n
```

This is the confidence sequence Splitcheck reports. It's simultaneously
valid at every sample size — you can peek every day and it still covers
the true lift with probability ≥ 1 - α.

**Choice of τ.** Bigger τ (bigger prior on effect size) gives more power at
large effects and less at small ones; smaller τ does the opposite. Default:
`τ = 0.02` (2 percentage points absolute) for the analyzer, `τ = 0.01` for
the peeking sim. There is no "right" τ — it's a modeling choice, like the
MDE in a sample-size calculation.

**Peek-tolerance tax.** The always-valid CI is wider than the fixed-horizon
CI, especially at small `n`. On the shipped sample data, it's about 60%
wider. That's the price of anytime-validity — no free lunch.

Implementation: [`stats.js:alwaysValidInference`](stats.js) /
[`analyze.py:always_valid`](analyze.py).

**How to pick τ in practice.** τ is a *modeling* choice, not a truth. The
always-valid guarantee holds for any τ > 0 — only power depends on it.

- **Match your MDE.** If you'd have planned the experiment around a 1
  percentage-point absolute lift, set `τ ≈ 0.01`. Effects near your MDE
  detect with something close to the fixed-horizon power.
- **Rule of thumb: `τ ≈ 5% · p̂_baseline`.** For a 5% baseline, `τ ≈ 0.0025`;
  for a 20% baseline, `τ ≈ 0.01`. Scales the prior with what "a meaningful
  lift" looks like at that baseline.
- **Too small.** The mixture concentrates near δ = 0, so genuinely large
  effects look surprising even under the alternative and Λ stays small.
  You lose power on the wins that matter most.
- **Too big.** The mixture spreads mass over effects you'd never see, and
  the small effects you do see get down-weighted. You lose power on realistic
  lifts.
- **Never set τ = 0.** That collapses the mixture to a point mass at H₀
  and Λ ≡ 1. The test never rejects.
- **Optimizely-style: pick τ per-metric from historical variance.** If
  you have prior experiments on the same metric, set τ to the standard
  deviation of past observed lifts. Splitcheck doesn't do this
  automatically — but it's what you'd do to tune it.

### 6. Sample Ratio Mismatch: chi-square

**The check.** For a target split `r_A : r_B` (say 0.5 : 0.5) and observed
counts `n_A, n_B` with total `N`, the expected counts are `E_A = N · r_A`,
`E_B = N · r_B`. Pearson chi-square:

```
                (n_A - E_A)²      (n_B - E_B)²
    χ²   =  ────────────  +  ────────────
                    E_A                E_B
```

Under the null "the assignment matches the target split," this is
approximately `χ²` distributed with `k - 1 = 1` degrees of freedom. Reject
if `χ² > χ²_{1, 1-p_thresh}`.

**Computing the tail without a special-function library.** For `df = 1`,
use the identity `χ²_1 = Z²` where `Z ~ N(0, 1)`:

```
P(χ²_1 > c)  =  P( |Z| > √c )  =  2 · (1 - Φ(√c))
```

Since Splitcheck already has a `normCdf`, this line is free. For `df > 1`,
Wilson–Hilferty: `(X/df)^(1/3)` is approximately normal with mean
`1 - 2/(9df)` and variance `2/(9df)`. Both are in the code.

**Kohavi's threshold.** Fire an alert at `p < 0.001`. At that level, SRM
almost certainly reflects a real bug in the randomizer, the traffic
splitter, or the exposure/conversion logging — not sampling noise. Do not
trust the metric analysis until you find and fix it.

Implementation: [`stats.js:srmCheck`, `chiSqSurvival1`](stats.js) /
[`analyze.py:srm_check`](analyze.py).

---

## Worked example: `sample_data.csv`

Everything above stays abstract until you see it on real numbers. Here's the
shipped 20,000-event synthetic dataset (baseline 5.0%, treatment 5.6%, seed
42) run through every formula by hand.

**Data.**

```
A: n_A = 10,000   x_A = 456   p̂_A = 0.04560
B: n_B = 10,000   x_B = 557   p̂_B = 0.05570
Δ̂ = p̂_B - p̂_A = 0.01010    (1.01 percentage points, 22.15% relative lift)
```

**SRM check.** Expected 50/50 split, so `E_A = E_B = 10,000` exactly.

```
χ² = (10000 - 10000)²/10000 + (10000 - 10000)²/10000  =  0
p  = 2·(1 - Φ(0))  =  1
```

Assignment OK. Proceed.

**Frequentist z-test.** Pooled proportion, pooled SE for the test:

```
p̂    = (456 + 557) / 20000                          =  0.05065
SE₀  = √[ 0.05065 · 0.94935 · (1/10000 + 1/10000) ]  =  0.003101
z    = 0.01010 / 0.003101                            =  3.257
p    = 2·(1 - Φ(3.257))                              =  0.001126
```

Unpooled SE for the CI:

```
SE_Δ = √[ 0.04560·0.95440/10000 + 0.05570·0.94430/10000 ]  =  0.003100
CI   = 0.01010 ± 1.960 · 0.003100                            =  [0.402%, 1.618%]
```

Verdict at α = 0.05: **SIGNIFICANT.**

**Always-valid (mSPRT), τ = 0.02.**

```
V     = SE_Δ²  =  9.61 × 10⁻⁶
τ²    = 4.00 × 10⁻⁴
V/(V+τ²)     =  0.02347
log Λ = ½·log(0.02347) + (0.01010² · τ²) / (2 V (V+τ²))
      = -1.876 + 5.183
      =  3.307
p*    = e^(-3.307)  =  0.0367
```

Confidence sequence half-width:

```
w = √[ 2V(V+τ²)/τ² · ( log(1/0.05) + ½·log((V+τ²)/V) ) ]
  = √[ 1.968 × 10⁻⁵ · (2.996 + 1.876) ]
  = 0.00979
CS = 0.01010 ± 0.00979  =  [0.031%, 1.989%]
```

The CS is ~60% wider than the fixed-horizon CI — that's the peek-tolerance tax.

**Bayesian.** Posteriors (flat prior):

```
p_A | data ~ Beta(1 + 456, 1 + 9544)  =  Beta(457, 9545)
p_B | data ~ Beta(1 + 557, 1 + 9443)  =  Beta(558, 9444)
```

Monte Carlo (20k draws):

```
P(B > A)                =  100.0%     (all draws had p_B > p_A)
E[loss | pick A]        =  1.011%     ≈ the true lift
E[loss | pick B]        =  0.000%
95% CrI on p_A          =  [4.17%, 4.99%]
95% CrI on p_B          =  [5.14%, 6.04%]
```

All four frames agree: **B wins.** Numbers here match the Python CLI
(`python3 analyze.py sample_data.csv --seed 1`) to the last printed digit.

---

## Design decisions

- **No dependencies.** `stats.js` is one file, `analyze.py` is one file.
  Both are 200-ish lines. Everything numeric is written from scratch so
  the derivations line up with the code.
- **Split-standard-error convention** for the frequentist test (pooled for
  p-value, unpooled for CI). This is what R's `prop.test` and Kohavi's
  book use. Common tutorial mistake avoided.
- **Beta(1, 1) prior** for the Bayesian analysis. Uniform on [0, 1], no
  strong prior belief. If you have historical data, plug in different
  hyperparameters.
- **Gaussian mixing prior** for the mSPRT, not a point alternative or a
  discrete mixture. Analytic, one hyperparameter (τ), works for any effect
  size range.
- **Ville's inequality as the peek-tolerance guarantee**, not alpha-spending
  or O'Brien-Fleming boundaries. Cleaner theory, simpler code, doesn't
  require you to commit to interim-look times.
- **SRM front and center.** Every result view leads with it. If it fires,
  the metric analysis is untrustworthy and there's no point looking at
  the rest.

## What Splitcheck is NOT

Being honest about the scope:

- **Not an experimentation platform.** No user bucketing, no feature-flag
  delivery, no event pipeline, no dashboards. If you want Statsig /
  GrowthBook / Optimizely, this is not that.
- **Not multi-armed.** Two arms only. No multi-armed bandits, no factorial
  designs, no interaction analysis.
- **Not multi-metric.** One conversion metric per experiment. No Bonferroni
  or Benjamini–Hochberg across families of metrics.
- **Not CUPED.** No variance reduction from pre-experiment covariates.
- **Not for continuous metrics.** Two-proportion machinery. For revenue,
  session duration, etc. you'd need a t-test or a bootstrap.
- **Not for very small samples.** The normal approximation assumes each
  cell has at least ~5–10 conversions. For rare-event metrics use Fisher
  exact or a Beta-Binomial exact calculation.

## Files

| File | Purpose |
| --- | --- |
| `index.html` | Three-tab UI. |
| `app.js` | UI wiring — reads inputs, calls `stats.js`, renders results. |
| `stats.js` | All statistics. No dependencies. |
| `styles.css` | Dark theme. |
| `analyze.py` | Python CLI. Reads `variant,converted` CSV, prints all four analyses. |
| `sample_data.csv` | 20,000 synthetic events (A: 5.0%, B: 5.6%). Seed 42. |

## References

**Two-proportion z-test and sample size:**
- Wikipedia, [Two-proportion Z-test](https://en.wikipedia.org/wiki/Two-proportion_Z-test).
- Fleiss, Levin, Paik. *Statistical Methods for Rates and Proportions*, 3rd ed. Wiley.

**Bayesian A/B testing:**
- Chris Stucchio (2015). [Bayesian A/B Testing at VWO](https://cdn2.hubspot.net/hubfs/310840/VWO_SmartStats_technical_whitepaper.pdf).
- Michael Frasco. [The Power of Bayesian A/B Testing](https://medium.com/convoy-tech/the-power-of-bayesian-a-b-testing-f859d2219d5), Convoy Tech.
- David Robinson (2015). [Is Bayesian A/B Testing Immune to Peeking? Not Exactly](http://varianceexplained.org/r/bayesian-ab-testing/).

**Beta / Gamma sampling:**
- Marsaglia & Tsang (2000). "A Simple Method for Generating Gamma Variables." *ACM TOMS* 26(3).

**Always-valid inference / mSPRT:**
- Johari, Koomen, Pekelis, Walsh (2015/2022). "Always Valid Inference: Continuous Monitoring of A/B Tests." *Operations Research*.
- Ramesh Johari et al. (2017). ["Peeking at A/B Tests." KDD](http://library.usc.edu.ph/ACM/KKD%202017/pdfs/p1517.pdf).
- Howard, Ramdas, McAuliffe, Sekhon. ["Time-uniform, nonparametric, nonasymptotic confidence sequences." *Annals of Statistics* (2021)](https://arxiv.org/abs/1810.08240).
- Ville (1939). *Étude critique de la notion de collectif.* Original martingale-inequality paper.

**Production experimentation methodology:**
- Kohavi, Tang, Xu. *Trustworthy Online Controlled Experiments.* Cambridge University Press, 2020.
- Fabijan et al. (2019). ["Diagnosing Sample Ratio Mismatch in Online Controlled Experiments." KDD.](https://dl.acm.org/doi/10.1145/3292500.3330722)
- Spotify Engineering. [Choosing a Sequential Testing Framework](https://engineering.atspotify.com/2023/03/choosing-sequential-testing-framework-comparisons-and-discussions).
