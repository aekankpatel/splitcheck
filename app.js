import {
  sampleSizePerArm,
  analyzeFrequentist,
  analyzeBayesian,
  alwaysValidInference,
  srmCheck,
  simulatePeeking,
} from "./stats.js";

// ---- Tabs ------------------------------------------------------------------
document.querySelectorAll(".tab").forEach((btn) => {
  btn.addEventListener("click", () => {
    document
      .querySelectorAll(".tab")
      .forEach((b) => b.classList.remove("active"));
    document
      .querySelectorAll(".view")
      .forEach((v) => v.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById(btn.dataset.tab).classList.add("active");
  });
});

// ---- Formatting helpers ----------------------------------------------------
const pct = (x, d = 2) => (x * 100).toFixed(d) + "%";
const num = (x, d = 3) => {
  if (!Number.isFinite(x)) return "—";
  if (Math.abs(x) >= 1000) return x.toFixed(0);
  return x.toFixed(d);
};
const int = (x) => x.toLocaleString();

// ---- Sample size calculator ------------------------------------------------
document.getElementById("size-run").addEventListener("click", () => {
  const p1 = parseFloat(document.getElementById("size-baseline").value) / 100;
  const mde = parseFloat(document.getElementById("size-mde").value) / 100;
  const alpha = parseFloat(document.getElementById("size-alpha").value);
  const power = parseFloat(document.getElementById("size-power").value);
  const out = document.getElementById("size-out");
  try {
    const n = sampleSizePerArm({ p1, mde, alpha, power });
    const abs = p1 * mde;
    out.classList.remove("hidden");
    out.innerHTML = `
      <div class="stat-row"><span class="k">Sample size per arm</span><span class="v"><strong>${int(n)}</strong></span></div>
      <div class="stat-row"><span class="k">Total sample size (both arms)</span><span class="v">${int(n * 2)}</span></div>
      <div class="stat-row"><span class="k">Baseline conversion (p1)</span><span class="v">${pct(p1)}</span></div>
      <div class="stat-row"><span class="k">Target conversion (p2)</span><span class="v">${pct(p1 * (1 + mde))}</span></div>
      <div class="stat-row"><span class="k">Absolute lift you can detect</span><span class="v">${pct(abs)}</span></div>
      <div class="explain">
        With <strong>${int(n)} visitors per arm</strong>, if the true relative lift is
        <strong>${pct(mde, 1)}</strong> or larger you'll detect it
        <strong>${pct(power, 0)}</strong> of the time. If there's really no effect, you'll
        <em>falsely</em> claim a winner only <strong>${pct(alpha, 0)}</strong> of the time.
        Smaller lifts need dramatically more traffic — halving the MDE roughly
        quadruples the required sample.
      </div>`;
  } catch (e) {
    out.classList.remove("hidden");
    out.innerHTML = `<div class="stat-row"><span class="k">Error</span><span class="v no">${e.message}</span></div>`;
  }
});

// ---- Frequentist + Bayesian analyzer ---------------------------------------
document.getElementById("analyze-run").addEventListener("click", () => {
  const nA = parseInt(document.getElementById("nA").value, 10);
  const xA = parseInt(document.getElementById("xA").value, 10);
  const nB = parseInt(document.getElementById("nB").value, 10);
  const xB = parseInt(document.getElementById("xB").value, 10);
  const alpha = parseFloat(document.getElementById("analyze-alpha").value);
  const splitA = parseFloat(document.getElementById("analyze-split").value) / 100;
  const tau = parseFloat(document.getElementById("analyze-tau").value);
  const out = document.getElementById("analyze-out");

  if ([nA, xA, nB, xB].some((v) => !Number.isFinite(v) || v < 0)) {
    out.classList.remove("hidden");
    out.innerHTML = `<div class="stat-row"><span class="k">Error</span><span class="v no">Enter non-negative integers.</span></div>`;
    return;
  }
  if (xA > nA || xB > nB) {
    out.classList.remove("hidden");
    out.innerHTML = `<div class="stat-row"><span class="k">Error</span><span class="v no">Conversions can't exceed visitors.</span></div>`;
    return;
  }

  const srm = srmCheck({ counts: [nA, nB], expected: [splitA, 1 - splitA] });
  const f = analyzeFrequentist({ nA, xA, nB, xB, alpha });
  const avi = alwaysValidInference({ nA, xA, nB, xB, tau, alpha });
  const bayes = analyzeBayesian({ nA, xA, nB, xB, samples: 20000 });

  const winner =
    f.pValue < alpha
      ? f.absLift > 0
        ? `<span class="badge ok">B wins</span>`
        : `<span class="badge ok">A wins</span>`
      : `<span class="badge neutral">Inconclusive</span>`;

  const srmBadge = srm.alert
    ? `<span class="badge no">SRM p=${srm.pValue < 1e-4 ? srm.pValue.toExponential(1) : srm.pValue.toFixed(4)} — don't trust results</span>`
    : `<span class="badge ok">Assignment OK (p=${srm.pValue.toFixed(3)})</span>`;

  out.classList.remove("hidden");
  out.innerHTML = `
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px;">
      <strong>Verdict</strong>
      ${winner}
    </div>

    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:12px;">
      <span style="font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:0.04em;">Sample Ratio Mismatch (χ²)</span>
      ${srmBadge}
    </div>
    ${
      srm.alert
        ? `<div class="explain" style="border-left-color: var(--danger); background: rgba(247,118,142,0.06); margin-top:0; margin-bottom:12px;">
             Observed split ${pct(srm.observed[0], 2)} / ${pct(srm.observed[1], 2)} deviates from expected
             ${pct(splitA, 1)} / ${pct(1 - splitA, 1)} more than random chance would allow
             (χ² = ${num(srm.chiSq, 2)}). Almost certainly a bug in the randomizer, bucketing, or logging.
             Fix that before trusting the metric analysis below — Kohavi's rule.
           </div>`
        : ""
    }

    <div class="two-col">
      <div class="col a">
        <h3>Variant A (control)</h3>
        <div class="stat-row"><span class="k">Visitors</span><span class="v">${int(nA)}</span></div>
        <div class="stat-row"><span class="k">Conversions</span><span class="v">${int(xA)}</span></div>
        <div class="stat-row"><span class="k">Rate</span><span class="v a">${pct(f.rateA, 3)}</span></div>
        <div class="stat-row"><span class="k">95% posterior on true rate</span><span class="v">${pct(bayes.posteriorA[0], 2)} – ${pct(bayes.posteriorA[2], 2)}</span></div>
      </div>
      <div class="col b">
        <h3>Variant B (treatment)</h3>
        <div class="stat-row"><span class="k">Visitors</span><span class="v">${int(nB)}</span></div>
        <div class="stat-row"><span class="k">Conversions</span><span class="v">${int(xB)}</span></div>
        <div class="stat-row"><span class="k">Rate</span><span class="v b">${pct(f.rateB, 3)}</span></div>
        <div class="stat-row"><span class="k">95% posterior on true rate</span><span class="v">${pct(bayes.posteriorB[0], 2)} – ${pct(bayes.posteriorB[2], 2)}</span></div>
      </div>
    </div>

    <h3 style="font-size:13px;color:var(--muted);text-transform:uppercase;margin:16px 0 4px;">Frequentist z-test</h3>
    <div class="stat-row"><span class="k">Absolute lift (B - A)</span><span class="v">${pct(f.absLift, 3)}</span></div>
    <div class="stat-row"><span class="k">Relative lift</span><span class="v">${Number.isFinite(f.relLift) ? pct(f.relLift, 2) : "—"}</span></div>
    <div class="stat-row"><span class="k">z-score</span><span class="v">${num(f.z, 3)}</span></div>
    <div class="stat-row"><span class="k">p-value (two-sided)</span><span class="v ${f.pValue < alpha ? "ok" : "no"}">${f.pValue < 1e-4 ? f.pValue.toExponential(2) : f.pValue.toFixed(4)}</span></div>
    <div class="stat-row"><span class="k">95% CI on lift</span><span class="v">[${pct(f.ci95[0], 3)}, ${pct(f.ci95[1], 3)}]</span></div>

    <h3 style="font-size:13px;color:var(--muted);text-transform:uppercase;margin:16px 0 4px;">Always-valid (mSPRT, Johari et al.)</h3>
    <div class="stat-row"><span class="k">Always-valid p-value</span><span class="v ${avi.pValue < alpha ? "ok" : "no"}">${avi.pValue < 1e-4 ? avi.pValue.toExponential(2) : avi.pValue.toFixed(4)}</span></div>
    <div class="stat-row"><span class="k">Confidence sequence on lift</span><span class="v">[${pct(avi.ci[0], 3)}, ${pct(avi.ci[1], 3)}]</span></div>
    <div class="stat-row"><span class="k">Prior scale τ</span><span class="v">${pct(tau, 2)} absolute</span></div>
    <div class="stat-row"><span class="k">Extra width vs fixed-horizon CI</span><span class="v">${(((avi.ci[1] - avi.ci[0]) / (f.ci95[1] - f.ci95[0])) * 100 - 100).toFixed(0)}%</span></div>

    <h3 style="font-size:13px;color:var(--muted);text-transform:uppercase;margin:16px 0 4px;">Bayesian (Beta-Binomial, flat prior)</h3>
    <div class="stat-row"><span class="k">P(B beats A)</span><span class="v ${bayes.probBBeatsA > 0.95 || bayes.probBBeatsA < 0.05 ? "ok" : ""}">${pct(bayes.probBBeatsA, 1)}</span></div>
    <div class="stat-row"><span class="k">Expected loss if you pick A</span><span class="v">${pct(bayes.expectedLossPickingA, 3)}</span></div>
    <div class="stat-row"><span class="k">Expected loss if you pick B</span><span class="v">${pct(bayes.expectedLossPickingB, 3)}</span></div>

    <div class="explain">
      The <strong>fixed-horizon p-value</strong> is only valid if you decided the sample size
      before looking. The <strong>always-valid p-value</strong> is safe under continuous monitoring
      — you can peek every day and still control the false-positive rate at α. It's wider
      early on and catches up as the experiment matures. The Bayesian <strong>P(B beats A)</strong>
      is a direct posterior probability, and <strong>expected loss</strong> tells you what
      conversion you'd give up if you picked the wrong arm. <strong>SRM</strong> at the top is
      the trust check — if it fires, don't believe any of the numbers below it until you fix
      the assignment bug.
    </div>
  `;
});

// ---- Peeking simulator -----------------------------------------------------
document.getElementById("peek-run").addEventListener("click", async () => {
  const rate = parseFloat(document.getElementById("peek-rate").value) / 100;
  const nPerArm = parseInt(document.getElementById("peek-n").value, 10);
  const checkEvery = parseInt(document.getElementById("peek-freq").value, 10);
  const trials = parseInt(document.getElementById("peek-trials").value, 10);
  const alpha = parseFloat(document.getElementById("peek-alpha").value);
  const tau = parseFloat(document.getElementById("peek-tau").value);
  const out = document.getElementById("peek-out");
  const btn = document.getElementById("peek-run");

  btn.disabled = true;
  btn.textContent = "Simulating…";
  out.classList.remove("hidden");
  out.innerHTML = `<div class="stat-row"><span class="k">Running ${int(trials)} trials…</span><span class="v">—</span></div>`;

  // Yield so the button state can paint before the heavy loop runs.
  await new Promise((r) => setTimeout(r, 20));

  const res = simulatePeeking({
    trueRate: rate,
    nPerArm,
    checkEvery,
    alpha,
    trials,
    tau,
  });

  const inflation =
    res.honestFalsePositiveRate === 0
      ? Infinity
      : res.peekingFalsePositiveRate / res.honestFalsePositiveRate;

  out.innerHTML = `
    <div class="stat-row"><span class="k">Honest fixed-horizon (one test at the end)</span>
      <span class="v ok">${pct(res.honestFalsePositiveRate, 1)}</span></div>
    <div class="stat-row"><span class="k">Naive peeker (fixed-horizon, stop at first p&lt;${alpha})</span>
      <span class="v no">${pct(res.peekingFalsePositiveRate, 1)}</span></div>
    <div class="stat-row"><span class="k">Always-valid peeker (mSPRT, τ=${pct(tau, 2)})</span>
      <span class="v ok">${pct(res.alwaysValidPeekerFalsePositiveRate, 1)}</span></div>
    <div class="stat-row"><span class="k">Naive-peeker inflation vs α</span>
      <span class="v no">${Number.isFinite(inflation) ? inflation.toFixed(1) + "×" : "∞"}</span></div>

    <div class="explain">
      Both arms have the same true rate, so any "winner" is a false positive. The <strong>honest</strong>
      test at n=${int(nPerArm)} hits the nominal <strong>${pct(alpha, 0)}</strong> — expected.
      The <strong>naive peeker</strong> checks every ${int(checkEvery)} visitors and stops on the first
      p&lt;${alpha}; with ${int(nPerArm / checkEvery)} looks, random noise crosses the line far more
      than ${pct(alpha, 0)} of the time — that's the peeking problem. The
      <strong>always-valid peeker</strong> peeks just as often but uses the mSPRT p-value (a
      martingale under H₀), so Ville's inequality bounds the false-positive rate at α no matter
      how many times you look. This is the fix — you pay for it in slightly wider CIs.
    </div>
  `;

  btn.disabled = false;
  btn.textContent = "Run simulation";
});
