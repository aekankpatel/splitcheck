// Statistics for A/B testing (two-proportion z-test + Beta-Binomial Bayesian).
// No external deps. All functions are pure.

// Abramowitz & Stegun 26.2.17 — normal CDF, good to ~7.5e-8.
export function normCdf(z) {
  const sign = z < 0 ? -1 : 1;
  const x = Math.abs(z) / Math.SQRT2;
  const t = 1 / (1 + 0.3275911 * x);
  const y =
    1 -
    ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) *
      t +
      0.254829592) *
      t *
      Math.exp(-x * x);
  return 0.5 * (1 + sign * y);
}

// Inverse normal CDF via Beasley-Springer-Moro. Accurate enough for z_alpha/z_beta.
export function normInv(p) {
  if (p <= 0 || p >= 1) throw new Error("normInv: p must be in (0,1)");
  const a = [
    -39.6968302866538, 220.946098424521, -275.928510446969, 138.357751867269,
    -30.6647980661472, 2.50662827745924,
  ];
  const b = [
    -54.4760987982241, 161.585836858041, -155.698979859887, 66.8013118877197,
    -13.2806815528857,
  ];
  const c = [
    -0.00778489400243029, -0.322396458041136, -2.40075827716184,
    -2.54973253934373, 4.37466414146497, 2.93816398269878,
  ];
  const d = [
    0.00778469570904146, 0.32246712907004, 2.445134137143, 3.75440866190742,
  ];
  const plow = 0.02425;
  const phigh = 1 - plow;
  let q, r;
  if (p < plow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (
      (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    );
  }
  if (p <= phigh) {
    q = p - 0.5;
    r = q * q;
    return (
      ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) *
        q) /
      (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
    );
  }
  q = Math.sqrt(-2 * Math.log(1 - p));
  return (
    -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
    ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
  );
}

// Chi-square df=1 survival function via χ²_1 = Z², so
// P(χ²_1 > c) = P(|Z| > √c) = 2·(1 - Φ(√c)).
export function chiSqSurvival1(x) {
  if (x <= 0) return 1;
  return 2 * (1 - normCdf(Math.sqrt(x)));
}

// Wilson-Hilferty normal approximation for higher df (unused for 2-arm SRM, but
// there in case you add multi-arm later).
function chiSqSurvival(x, df) {
  if (df === 1) return chiSqSurvival1(x);
  if (x <= 0) return 1;
  const t = Math.pow(x / df, 1 / 3);
  const mu = 1 - 2 / (9 * df);
  const sigma = Math.sqrt(2 / (9 * df));
  return 1 - normCdf((t - mu) / sigma);
}

// Sample Ratio Mismatch check. counts is [nA, nB, ...] observed; expected is
// [rA, rB, ...] target ratios that sum to 1. Kohavi's rule of thumb: alert at
// p < 0.001, since SRM at that level almost certainly means the randomizer,
// bucketing, or logging has a bug and downstream metric results can't be trusted.
export function srmCheck({ counts, expected }) {
  const total = counts.reduce((a, b) => a + b, 0);
  if (total === 0) {
    return {
      chiSq: 0,
      pValue: 1,
      alert: false,
      total: 0,
      observed: [],
      expectedCounts: [],
    };
  }
  let chi = 0;
  const expectedCounts = [];
  const observed = [];
  for (let i = 0; i < counts.length; i++) {
    const e = total * expected[i];
    expectedCounts.push(e);
    observed.push(counts[i] / total);
    if (e > 0) chi += (counts[i] - e) ** 2 / e;
  }
  const df = counts.length - 1;
  return {
    chiSq: chi,
    pValue: chiSqSurvival(chi, df),
    alert: chiSqSurvival(chi, df) < 0.001,
    total,
    observed,
    expectedCounts,
  };
}

// Required sample size per arm for two-proportion z-test.
// p1: baseline conversion, mde: relative lift (e.g. 0.05 = 5% relative lift).
// alpha: two-sided sig level, power: 1 - beta.
export function sampleSizePerArm({ p1, mde, alpha = 0.05, power = 0.8 }) {
  const p2 = p1 * (1 + mde);
  if (p2 <= 0 || p2 >= 1) throw new Error("p2 out of range — reduce MDE");
  const zA = normInv(1 - alpha / 2);
  const zB = normInv(power);
  const num = Math.pow(zA + zB, 2) * (p1 * (1 - p1) + p2 * (1 - p2));
  const den = Math.pow(p2 - p1, 2);
  return Math.ceil(num / den);
}

// Two-proportion z-test on completed experiment.
// Returns {rateA, rateB, absLift, relLift, z, pValue, ci95, significant}.
export function analyzeFrequentist({ nA, xA, nB, xB, alpha = 0.05 }) {
  if (nA <= 0 || nB <= 0) throw new Error("sample sizes must be > 0");
  const pA = xA / nA;
  const pB = xB / nB;
  const pPool = (xA + xB) / (nA + nB);
  const sePool = Math.sqrt(pPool * (1 - pPool) * (1 / nA + 1 / nB));
  const z = sePool === 0 ? 0 : (pB - pA) / sePool;
  const pValue = 2 * (1 - normCdf(Math.abs(z)));
  const seDiff = Math.sqrt((pA * (1 - pA)) / nA + (pB * (1 - pB)) / nB);
  const zA = normInv(1 - alpha / 2);
  const ci95 = [pB - pA - zA * seDiff, pB - pA + zA * seDiff];
  return {
    rateA: pA,
    rateB: pB,
    absLift: pB - pA,
    relLift: pA === 0 ? NaN : (pB - pA) / pA,
    z,
    pValue,
    ci95,
    significant: pValue < alpha,
  };
}

// Always-valid inference via mSPRT with a Gaussian mixing prior N(0, τ²) on the
// true difference in proportions (Johari, Koomen, Pekelis, Walsh 2015, "Always
// Valid Inference"). Gives a p-value and confidence interval that are simultaneously
// valid at every sample size — so you can peek as often as you like without
// inflating the false-positive rate. Cost: intervals are wider than the fixed-
// horizon ones, especially early on.
//
// tau is your prior scale for the effect size (absolute lift). Bigger tau → more
// power at big effects, less at small ones. Default 0.02 = 2 percentage points.
export function alwaysValidInference({ nA, xA, nB, xB, tau = 0.02, alpha = 0.05 }) {
  if (nA <= 0 || nB <= 0) throw new Error("sample sizes must be > 0");
  const pA = xA / nA;
  const pB = xB / nB;
  const delta = pB - pA;
  const V = (pA * (1 - pA)) / nA + (pB * (1 - pB)) / nB;
  if (V === 0) {
    return { pValue: 1, ci: [delta, delta], halfWidth: 0, tau };
  }
  const t2 = tau * tau;
  // Log-likelihood ratio under the Gaussian mixture; testing H0: δ = 0.
  const logLambda =
    0.5 * Math.log(V / (V + t2)) +
    (delta * delta * t2) / (2 * V * (V + t2));
  const pValue = Math.min(1, Math.exp(-logLambda));
  // Confidence sequence: invert the mSPRT at level alpha for every candidate δ.
  const halfWidth = Math.sqrt(
    ((2 * V * (V + t2)) / t2) *
      (Math.log(1 / alpha) + 0.5 * Math.log((V + t2) / V))
  );
  return {
    pValue,
    ci: [delta - halfWidth, delta + halfWidth],
    halfWidth,
    tau,
  };
}

// Sample from Beta(alpha, beta) via two Gammas — Marsaglia-Tsang.
function sampleGamma(shape) {
  if (shape < 1) {
    const u = Math.random();
    return sampleGamma(shape + 1) * Math.pow(u, 1 / shape);
  }
  const d = shape - 1 / 3;
  const c = 1 / Math.sqrt(9 * d);
  while (true) {
    let x, v;
    do {
      // Box-Muller for a standard normal
      const u1 = Math.random();
      const u2 = Math.random();
      x = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
      v = 1 + c * x;
    } while (v <= 0);
    v = v * v * v;
    const u = Math.random();
    if (u < 1 - 0.0331 * x * x * x * x) return d * v;
    if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
  }
}
function sampleBeta(a, b) {
  const x = sampleGamma(a);
  const y = sampleGamma(b);
  return x / (x + y);
}

// Bayesian analysis with Beta(1,1) prior. P(B > A), expected loss, credible intervals.
export function analyzeBayesian({ nA, xA, nB, xB, samples = 20000 }) {
  const aA = 1 + xA;
  const bA = 1 + nA - xA;
  const aB = 1 + xB;
  const bB = 1 + nB - xB;
  let wins = 0;
  let sumLossA = 0; // expected loss of picking A when B is truly better
  let sumLossB = 0;
  const drawsA = new Float64Array(samples);
  const drawsB = new Float64Array(samples);
  for (let i = 0; i < samples; i++) {
    const a = sampleBeta(aA, bA);
    const b = sampleBeta(aB, bB);
    drawsA[i] = a;
    drawsB[i] = b;
    if (b > a) {
      wins++;
      sumLossA += b - a;
    } else {
      sumLossB += a - b;
    }
  }
  return {
    probBBeatsA: wins / samples,
    expectedLossPickingA: sumLossA / samples,
    expectedLossPickingB: sumLossB / samples,
    posteriorA: quantiles(drawsA, [0.025, 0.5, 0.975]),
    posteriorB: quantiles(drawsB, [0.025, 0.5, 0.975]),
  };
}

function quantiles(arr, qs) {
  const sorted = Array.from(arr).sort((a, b) => a - b);
  return qs.map((q) => sorted[Math.floor(q * (sorted.length - 1))]);
}

// Peeking-problem simulator. Under H0 (both arms identical), compares three
// strategies: an honest experimenter (one test at nPerArm), a naive peeker
// (repeated fixed-horizon tests, stop on first p<alpha), and an always-valid
// peeker (repeated mSPRT tests). The first should hit ~alpha; the second
// inflates badly; the third should also hit ~alpha because always-valid
// inference is peek-tolerant by construction.
export function simulatePeeking({
  trueRate = 0.1,
  nPerArm = 5000,
  checkEvery = 100,
  alpha = 0.05,
  trials = 500,
  tau = 0.01,
}) {
  let naivePeekerFP = 0;
  let alwaysValidPeekerFP = 0;
  let honestFP = 0;
  for (let t = 0; t < trials; t++) {
    let xA = 0;
    let xB = 0;
    let naiveFlagged = false;
    let aviFlagged = false;
    for (let i = 1; i <= nPerArm; i++) {
      if (Math.random() < trueRate) xA++;
      if (Math.random() < trueRate) xB++;
      if (i % checkEvery === 0) {
        if (!naiveFlagged) {
          const { pValue } = analyzeFrequentist({ nA: i, xA, nB: i, xB, alpha });
          if (pValue < alpha) naiveFlagged = true;
        }
        if (!aviFlagged) {
          const { pValue } = alwaysValidInference({ nA: i, xA, nB: i, xB, tau, alpha });
          if (pValue < alpha) aviFlagged = true;
        }
      }
    }
    if (naiveFlagged) naivePeekerFP++;
    if (aviFlagged) alwaysValidPeekerFP++;
    const final = analyzeFrequentist({ nA: nPerArm, xA, nB: nPerArm, xB, alpha });
    if (final.pValue < alpha) honestFP++;
  }
  return {
    peekingFalsePositiveRate: naivePeekerFP / trials,
    alwaysValidPeekerFalsePositiveRate: alwaysValidPeekerFP / trials,
    honestFalsePositiveRate: honestFP / trials,
    trials,
  };
}
