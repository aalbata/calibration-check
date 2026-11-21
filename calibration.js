/*
 * Calibration measures, calibration maps, delayed refit schedules and the decision
 * procedure of the article, written to follow temporal_calibration/ (Python) step by
 * step. Runs in a browser (global TC) and in Node (module.exports). No dependencies.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) {
    module.exports = factory();
  } else {
    root.TC = factory();
  }
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  var PCLIP = 1e-6;
  var SPARSE_LEFT_OUT = 0.05;
  var ROLE_ORDER = { calibration: 0, development: 1, deployment: 2 };
  var FAMILIES = ["uncalibrated", "platt", "temperature"];
  var FAMILY_LABEL = {
    uncalibrated: "uncalibrated",
    platt: "Platt scaling",
    temperature: "temperature scaling"
  };

  function clip(p) {
    return p < PCLIP ? PCLIP : p > 1 - PCLIP ? 1 - PCLIP : p;
  }

  function expit(z) {
    if (z >= 0) {
      return 1 / (1 + Math.exp(-z));
    }
    var e = Math.exp(z);
    return e / (1 + e);
  }

  function logit(p) {
    return Math.log(p / (1 - p));
  }

  function mean(a) {
    var s = 0;
    for (var i = 0; i < a.length; i++) s += a[i];
    return a.length ? s / a.length : NaN;
  }

  /* Quantile with linear interpolation, as numpy.quantile (default method). */
  function quantile(sorted, q) {
    var n = sorted.length;
    var vi = n * q + (1 + q * (1 - 1 - 1)) - 1;
    var lo = Math.floor(vi);
    var g = vi - lo;
    if (lo < 0) { lo = 0; g = 0; }
    if (lo >= n - 1) return sorted[n - 1];
    var a = sorted[lo], b = sorted[lo + 1], d = b - a;
    return g >= 0.5 ? b - d * (1 - g) : a + d * g;
  }

  /* Number of edges strictly below p (numpy.searchsorted, side "left"). */
  function binIndex(edges, p) {
    var lo = 0, hi = edges.length;
    while (lo < hi) {
      var mid = (lo + hi) >> 1;
      if (edges[mid] < p) lo = mid + 1; else hi = mid;
    }
    return lo;
  }

  /* Equal-mass bins fixed from the fit-window outputs of a map; tied edges are merged. */
  function makeBinning(fitProbs, nBins) {
    nBins = nBins || 15;
    var sorted = Float64Array.from(fitProbs).sort();
    var step = 1 / nBins, edges = [];
    for (var i = 1; i < nBins; i++) {
      var e = quantile(sorted, i * step);
      if (!edges.length || e !== edges[edges.length - 1]) edges.push(e);
    }
    var b = edges.length + 1;
    var sum = new Float64Array(b), cnt = new Float64Array(b);
    for (var k = 0; k < fitProbs.length; k++) {
      var j = binIndex(edges, fitProbs[k]);
      sum[j] += fitProbs[k];
      cnt[j] += 1;
    }
    var values = new Float64Array(b);
    for (var m = 0; m < b; m++) {
      values[m] = cnt[m] > 0 ? sum[m] / cnt[m]
        : ((m > 0 ? edges[m - 1] : 0) + (m < b - 1 ? edges[m] : 1)) / 2;
    }
    return { edges: edges, values: values, nBins: b };
  }

  /* Debiased squared calibration error of the discretized map on fixed bins. */
  function dsce(probs, y, binning) {
    var b = binning.nBins, n = y.length;
    var cnt = new Float64Array(b), pos = new Float64Array(b);
    for (var i = 0; i < n; i++) {
      var j = binIndex(binning.edges, probs[i]);
      cnt[j] += 1;
      pos[j] += y[i];
    }
    var plugin = 0, deb = 0, degenerate = 0, nUsed = 0;
    for (var k = 0; k < b; k++) {
      if (cnt[k] <= 1) { degenerate += 1; continue; }
      var ph = cnt[k] / n, yb = pos[k] / cnt[k];
      var gap2 = (binning.values[k] - yb) * (binning.values[k] - yb);
      plugin += ph * gap2;
      deb += ph * (gap2 - yb * (1 - yb) / (cnt[k] - 1));
      nUsed += cnt[k];
    }
    return {
      dsce: Math.max(deb, 0), dsceUnfloored: deb, rmsce: Math.sqrt(Math.max(deb, 0)),
      plugin: plugin, floored: deb < 0, nBinsDegenerate: degenerate, n: n,
      nBins: b, nBinsUsed: b - degenerate, nUsed: nUsed
    };
  }

  /*
   * Whether a period holds enough samples for its calibration error to be read. A bin
   * with fewer than two samples does not enter the estimate, so a period can show an
   * error of zero only because its samples were left out. estimable is false when no
   * bin entered; sparse is true when the period holds fewer than two samples per bin
   * or when more than SPARSE_LEFT_OUT of its samples were left out.
   */
  function evidence(d) {
    var estimable = d.nBinsUsed > 0;
    var sparse = !estimable || d.n < 2 * d.nBins || d.nUsed < (1 - SPARSE_LEFT_OUT) * d.n;
    return { estimable: estimable, sparse: sparse };
  }

  /* Area under the ROC curve with tied probabilities at their midrank. */
  function auc(probs, y) {
    var n = y.length, idx = new Uint32Array(n);
    for (var i = 0; i < n; i++) idx[i] = i;
    idx.sort(function (a, b) { return probs[a] - probs[b]; });
    var nPos = 0, rankSum = 0, i0 = 0;
    while (i0 < n) {
      var i1 = i0;
      while (i1 + 1 < n && probs[idx[i1 + 1]] === probs[idx[i0]]) i1++;
      var mid = (i0 + i1) / 2 + 1;
      for (var t = i0; t <= i1; t++) {
        if (y[idx[t]] === 1) { nPos += 1; rankSum += mid; }
      }
      i0 = i1 + 1;
    }
    var nNeg = n - nPos;
    if (nPos === 0 || nNeg === 0) return NaN;
    return (rankSum - nPos * (nPos + 1) / 2) / (nPos * nNeg);
  }

  function prevalence(y) {
    var s = 0;
    for (var i = 0; i < y.length; i++) s += y[i];
    return y.length ? s / y.length : NaN;
  }

  function hasBothClasses(y) {
    var s = 0;
    for (var i = 0; i < y.length; i++) s += y[i];
    return s > 0 && s < y.length;
  }

  /*
   * Logistic regression of y on one feature x with intercept, minimizing the summed
   * log loss plus slope^2 / (2 C), the objective of scikit-learn with C = 1e6.
   * Newton iterations with step halving.
   */
  function fitLogistic(x, y, C) {
    C = C || 1e6;
    var n = y.length, a = 0, b = 0, lam = 1 / C;
    function loss(aa, bb) {
      var s = 0;
      for (var i = 0; i < n; i++) {
        var t = aa * x[i] + bb;
        s += (t > 0 ? t + Math.log1p(Math.exp(-t)) : Math.log1p(Math.exp(t))) - y[i] * t;
      }
      return s + 0.5 * lam * aa * aa;
    }
    var cur = loss(a, b);
    for (var it = 0; it < 200; it++) {
      var ga = lam * a, gb = 0, haa = lam, hab = 0, hbb = 0;
      for (var i = 0; i < n; i++) {
        var p = expit(a * x[i] + b), r = p - y[i], w = p * (1 - p);
        ga += r * x[i]; gb += r;
        haa += w * x[i] * x[i]; hab += w * x[i]; hbb += w;
      }
      var det = haa * hbb - hab * hab;
      if (!(det > 0)) break;
      var da = (hbb * ga - hab * gb) / det, db = (haa * gb - hab * ga) / det;
      var stepLen = 1, next = cur, na = a, nb = b;
      for (var h = 0; h < 40; h++) {
        na = a - stepLen * da; nb = b - stepLen * db;
        next = loss(na, nb);
        if (next <= cur) break;
        stepLen /= 2;
      }
      var moved = Math.abs(na - a) + Math.abs(nb - b);
      a = na; b = nb;
      if (moved < 1e-12 || Math.abs(cur - next) < 1e-13 * (1 + Math.abs(cur))) { cur = next; break; }
      cur = next;
    }
    return { slope: a, intercept: b };
  }

  /* Temperature T minimizing the log loss of expit(z / T), with log T in [-3, 3]. */
  function fitTemperature(z, y) {
    var n = y.length;
    function nll(logT) {
      var T = Math.exp(logT), s = 0;
      for (var i = 0; i < n; i++) {
        var p = clip(expit(z[i] / T));
        s -= y[i] ? Math.log(p) : Math.log(1 - p);
      }
      return s / n;
    }
    var gr = (Math.sqrt(5) - 1) / 2, lo = -3, hi = 3;
    var c = hi - gr * (hi - lo), d = lo + gr * (hi - lo), fc = nll(c), fd = nll(d);
    for (var it = 0; it < 200 && hi - lo > 1e-11; it++) {
      if (fc < fd) { hi = d; d = c; fd = fc; c = hi - gr * (hi - lo); fc = nll(c); }
      else { lo = c; c = d; fc = fd; d = lo + gr * (hi - lo); fd = nll(d); }
    }
    return Math.exp((lo + hi) / 2);
  }

  /* A fitted calibration map: predict(z) returns clipped probabilities. */
  function fitMap(name, z, y) {
    var map;
    if (name === "uncalibrated") {
      map = { name: name, params: {}, f: function (v) { return expit(v); } };
    } else if (name === "temperature") {
      var T = fitTemperature(z, y);
      map = { name: name, params: { T: T }, f: function (v) { return expit(v / T); } };
    } else if (name === "platt") {
      var lr = fitLogistic(z, y, 1e6);
      map = { name: name, params: { slope: lr.slope, intercept: lr.intercept },
              f: function (v) { return expit(lr.slope * v + lr.intercept); } };
    } else {
      throw new Error("unknown calibration map: " + name);
    }
    map.predict = function (zz) {
      var out = new Float64Array(zz.length);
      for (var i = 0; i < zz.length; i++) out[i] = clip(map.f(zz[i]));
      return out;
    };
    return map;
  }

  /* Calibration intercept and slope: logistic regression of y on logit(p). */
  function interceptSlope(probs, y) {
    if (!hasBothClasses(y)) return { intercept: NaN, slope: NaN };
    var x = new Float64Array(probs.length);
    for (var i = 0; i < probs.length; i++) x[i] = logit(clip(probs[i]));
    var lr = fitLogistic(x, y, 1e6);
    return { intercept: lr.intercept, slope: lr.slope };
  }

  /* Trapezoidal average over slots with equal slot weights. */
  function aut(v) {
    if (v.length === 1) return v[0];
    var s = 0;
    for (var i = 1; i < v.length; i++) s += (v[i] + v[i - 1]) / 2;
    return s / (v.length - 1);
  }

  /*
   * Fixed-cadence refits of one map family over ordered slots. At the start of slot k
   * the labels of slot j are usable when k >= j + 1 + delta; a refit happens when such
   * a slot exists and k is a multiple of the cadence, uses the latest usable slot, and
   * is skipped when that slot holds one class. cadence 0 keeps the first map.
   */
  function runCadence(slots, fitZ, fitY, name, cadence, delta) {
    var map0 = fitMap(name, fitZ, fitY);
    var binning = makeBinning(map0.predict(fitZ), 15);
    var map = map0, probs = [], refits = [], labels = 0, bought = {};
    for (var k = 0; k < slots.length; k++) {
      var j = k - 1 - delta;
      if (cadence > 0 && j >= 0 && k % cadence === 0) {
        if (!bought[j]) { bought[j] = true; labels += slots[j].y.length; }
        if (hasBothClasses(slots[j].y)) {
          map = fitMap(name, slots[j].z, slots[j].y);
          refits.push({ slot: k, fittedOn: j, params: map.params });
        }
      }
      probs.push(map.predict(slots[k].z));
    }
    return { probs: probs, binning: binning, refits: refits, labels: labels, map0: map0 };
  }

  function slotDsce(run, slots) {
    return slots.map(function (s, k) { return dsce(run.probs[k], s.y, run.binning); });
  }

  /*
   * Full analysis of one detector.
   * data = { calib: {z, y}, dev: [{label, z, y}], dep: [{label, z, y}] } with logit-scale
   * scores z and labels y in {0, 1}; options = { delta, family, smallSpread }.
   */
  function analyze(data, options) {
    var delta = options.delta == null ? 1 : options.delta;
    var cz = data.calib.z, cy = data.calib.y, dev = data.dev, dep = data.dep;
    if (!hasBothClasses(cy)) throw new Error("The calibration window needs samples of both classes.");
    if (dev.length < 1) throw new Error("At least one development period is needed.");

    // Candidate maps, each scored on bins from its own fit-window outputs.
    var candidates = FAMILIES.map(function (name) {
      var run = runCadence(dev, cz, cy, name, 0, delta);
      var d = slotDsce(run, dev).map(function (r) { return r.dsce; });
      return { name: name, devAutDsce: aut(d), params: run.map0.params };
    });
    var selected = options.family && options.family !== "auto" ? options.family : null;
    if (!selected) {
      var best = Infinity;
      candidates.forEach(function (c) {
        if (c.devAutDsce < best) { best = c.devAutDsce; selected = c.name; }
      });
    }

    function perSlot(slots, run) {
      return slots.map(function (s, k) {
        var d = dsce(run.probs[k], s.y, run.binning), cs = interceptSlope(run.probs[k], s.y);
        var ev = evidence(d);
        return {
          label: s.label, n: s.y.length, prevalence: prevalence(s.y),
          auc: auc(run.probs[k], s.y), dsce: d.dsce, rmsce: d.rmsce, floored: d.floored,
          intercept: cs.intercept, slope: cs.slope,
          nBins: d.nBins, nBinsUsed: d.nBinsUsed, nUsed: d.nUsed,
          estimable: ev.estimable, sparse: ev.sparse
        };
      });
    }

    var devNever = runCadence(dev, cz, cy, selected, 0, delta);
    var devRows = perSlot(dev, devNever);
    var devAutNever = aut(devRows.map(function (r) { return r.dsce; }));
    // A map without parameters cannot be refitted, so refit schedules do not apply to it.
    var refittable = selected !== "uncalibrated";
    var backtest = [1, 2, 3].map(function (c) {
      var run = runCadence(dev, cz, cy, selected, c, delta);
      var d = slotDsce(run, dev);
      return { cadence: c, refits: refittable ? run.refits.length : 0,
               devAutDsce: aut(d.map(function (r) { return r.dsce; })),
               notEstimable: d.filter(function (r) { return r.nBinsUsed === 0; }).length };
    });
    var cStar = backtest[0];
    backtest.forEach(function (b) { if (b.devAutDsce < cStar.devAutDsce) cStar = b; });

    var depRows = [], depSchedules = [];
    if (dep.length) {
      var depNever = runCadence(dep, cz, cy, selected, 0, delta);
      depRows = perSlot(dep, depNever);
      [0, 1, 2, 3].forEach(function (c) {
        var run = c === 0 ? depNever : runCadence(dep, cz, cy, selected, c, delta);
        var d = slotDsce(run, dep);
        var applicable = c === 0 || refittable;
        var estimable = d.map(function (r) { return r.nBinsUsed > 0; });
        depSchedules.push({
          cadence: c, applicable: applicable,
          refits: applicable ? run.refits.length : 0, labels: applicable ? run.labels : 0,
          dsce: d.map(function (r) { return r.dsce; }),
          rmsce: d.map(function (r) { return r.rmsce; }),
          estimable: estimable,
          nEstimable: estimable.filter(Boolean).length,
          autDsce: aut(d.map(function (r) { return r.dsce; })),
          autRmsce: aut(d.map(function (r) { return r.rmsce; }))
        });
      });
    }

    var calibPrev = prevalence(cy);
    var devPrev = devRows.map(function (r) { return r.prevalence; });
    var result = {
      delta: delta, selected: selected, refittable: refittable, candidates: candidates,
      devEstimable: devRows.filter(function (r) { return r.estimable; }).length,
      calib: { n: cy.length, prevalence: calibPrev, params: devNever.map0.params,
               nBins: devNever.binning.nBins },
      dev: devRows, devAutDsce: devAutNever, devAutRmsce: aut(devRows.map(function (r) { return r.rmsce; })),
      backtest: backtest, cStar: cStar.cadence, cStarDevAutDsce: cStar.devAutDsce,
      dep: depRows, depSchedules: depSchedules
    };
    result.decision = decide(result, options);
    return result;
  }

  /*
   * Decision procedure of the article with point estimates. The first deployment
   * period is the first matured slot.
   */
  function decide(r, options) {
    var steps = [], first = r.dep.length ? r.dep[0] : null;
    if (!first) {
      return { steps: steps, action: "none",
               text: "No deployment period with labels was supplied, so the procedure cannot be applied." };
    }
    var needed = r.dev.concat([first]);
    var question1 = "Is the AUC of the first matured slot below its development range?";
    var oneClass = needed.filter(function (d) { return !isFinite(d.auc); });
    if (oneClass.length) {
      steps.push({ question: question1, answer: null, state: "insufficient",
                   detail: "The AUC cannot be computed for " + listLabels(oneClass) + ", which " +
                           (oneClass.length === 1 ? "holds" : "hold") + " samples of one class only" });
      return { steps: steps, action: "insufficient",
               text: "Insufficient data to evaluate this step. Each development period and the first " +
                     "deployment period need samples of both classes." };
    }
    var devAuc = r.dev.map(function (d) { return d.auc; });
    var devDsce = r.dev.map(function (d) { return d.dsce; });
    var minAuc = Math.min.apply(null, devAuc), maxDsce = Math.max.apply(null, devDsce);

    var q1 = first.auc < minAuc;
    steps.push({ question: question1,
                 answer: q1,
                 detail: "AUC of " + first.label + " = " + fmt(first.auc, 4) +
                         "; lowest development AUC = " + fmt(minAuc, 4) });
    if (q1) {
      return { steps: steps, action: "retrain", evaluated: false,
               text: "Retrain the scorer, then restart the procedure." };
    }

    var question2 = "Does a delayed refit lower DSCE on the development slots, or does the DSCE of the first matured slot lie above the development range?";
    var unread = needed.filter(function (d) { return d.estimable === false; });
    var unreadRefit = (r.backtest || []).some(function (b) { return b.notEstimable > 0; });
    if (unread.length || unreadRefit) {
      steps.push({ question: question2, answer: null, state: "insufficient",
                   detail: "No calibration error can be computed for " +
                           (unread.length ? listLabels(unread) : "a development period under a refitted map") +
                           ", because no bin holds two or more of its samples" });
      return { steps: steps, action: "insufficient",
               text: "Insufficient data to evaluate this step. Each development period and the first " +
                     "deployment period need two or more samples in at least one bin." };
    }
    var noMap = r.selected === "uncalibrated";
    var backtestFavors = !noMap && r.cStarDevAutDsce < r.devAutDsce, above = first.dsce > maxDsce;
    var q2 = backtestFavors || above;
    steps.push({ question: question2,
                 answer: q2,
                 detail: "Development slot-averaged DSCE = " + fmt(r.devAutDsce, 5) +
                         (noMap ? "; the values are used as given, so there is no map to refit"
                                : " without refits and " + fmt(r.cStarDevAutDsce, 5) + " with refits every " +
                                  (r.cStar === 1 ? "slot" : r.cStar + " slots") + " at a label delay of " + r.delta) +
                         "; DSCE of " + first.label + " = " +
                         fmt(first.dsce, 5) + "; highest development DSCE = " + fmt(maxDsce, 5) });
    if (!q2) {
      return { steps: steps, action: "none", evaluated: true,
               text: "No routine refresh. Compare each newly matured slot with the development range." };
    }

    var devPrev = r.dev.map(function (d) { return d.prevalence; });
    var loDev = Math.min.apply(null, devPrev), hiDev = Math.max.apply(null, devPrev);
    var gap = Math.abs(logit(first.prevalence) - logit(r.calib.prevalence));
    var devSpread = logit(hiDev) - logit(loDev);
    var q3 = gap > devSpread;
    steps.push({ question: "Does the log-odds gap between calibration-window and first-matured-slot prevalence exceed the development log-odds spread?",
                 answer: q3,
                 detail: "Prevalence " + fmt(r.calib.prevalence, 3) + " in the calibration window and " +
                         fmt(first.prevalence, 3) + " in " + first.label + ", a log-odds gap of " + fmt(gap, 2) +
                         "; development prevalence " + fmt(loDev, 3) + " to " + fmt(hiDev, 3) +
                         ", a log-odds spread of " + fmt(devSpread, 2) });
    if (!q3 && noMap) {
      return { steps: steps, action: "no-map",
               text: "The procedure calls for a refit of the selected map on the latest matured slot, but " +
                     "no map is fitted when the values are used as given. Select temperature or Platt " +
                     "scaling to evaluate refits." };
    }
    if (!q3) {
      return { steps: steps, action: "refit-selected", evaluated: true,
               text: "Refit the selected map (" + FAMILY_LABEL[r.selected] + ") on the latest matured slot." };
    }

    var all = devPrev.concat([first.prevalence]);
    var spread = logit(Math.max.apply(null, all)) - logit(Math.min.apply(null, all));
    var cut = options.smallSpread;
    if (!(cut > 0)) {
      steps.push({ question: "Is the prevalence log-odds spread over the development and matured slots small?",
                   answer: null,
                   detail: "Log-odds spread = " + fmt(spread, 2) + ". Enter a cut point for a small spread to continue; the article does not supply one" });
      return { steps: steps, action: "undetermined",
               text: "The last question needs a cut point for a small prevalence spread, which must be set for the deployment." };
    }
    var q4 = spread <= cut;
    steps.push({ question: "Is the prevalence log-odds spread over the development and matured slots small?",
                 answer: q4,
                 detail: "Log-odds spread = " + fmt(spread, 2) + "; cut point = " + fmt(cut, 2) });
    if (q4) {
      return { steps: steps, action: "refit-intercept", evaluated: false,
               text: "Refit the intercept on recent labels, or apply a prior-shift correction." };
    }
    return { steps: steps, action: "slope-only", evaluated: true,
             text: "Use slope-only (temperature) refits, and keep the uncalibrated baseline as a candidate map." };
  }

  function fmt(x, d) {
    return isFinite(x) ? x.toFixed(d) : "n/a";
  }

  /* "a", "a and b" or "a, b and c" from the labels of periods. */
  function listLabels(rows) {
    var names = rows.map(function (d) { return d.label; });
    if (names.length < 3) return names.join(" and ");
    return names.slice(0, -1).join(", ") + " and " + names[names.length - 1];
  }

  var NUMBER = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;

  /*
   * Minimal CSV reader: comma, semicolon or tab separated, optional double quotes.
   * A row is kept only when its period is not empty, its label is 0 or 1 and its score
   * or probability is a written number; an empty cell is never read as zero. Rows left
   * out are counted in skipped and described in invalid with their line in the file.
   */
  function parseCsv(text) {
    var numbered = [];
    text.replace(/\r\n?/g, "\n").split("\n").forEach(function (l, i) {
      if (l.trim() !== "") numbered.push({ line: i + 1, text: l });
    });
    var lines = numbered.map(function (l) { return l.text; });
    if (lines.length < 2) throw new Error("The file needs a header row and at least one data row.");
    var sep = lines[0].indexOf("\t") >= 0 ? "\t" : lines[0].indexOf(";") >= 0 && lines[0].indexOf(",") < 0 ? ";" : ",";
    function split(line) {
      var out = [], cur = "", q = false;
      for (var i = 0; i < line.length; i++) {
        var ch = line[i];
        if (q) {
          if (ch === '"' && line[i + 1] === '"') { cur += '"'; i++; }
          else if (ch === '"') q = false;
          else cur += ch;
        } else if (ch === '"') q = true;
        else if (ch === sep) { out.push(cur); cur = ""; }
        else cur += ch;
      }
      out.push(cur);
      return out.map(function (s) { return s.trim(); });
    }
    var header = split(lines[0]).map(function (h) { return h.toLowerCase(); });
    function col(names) {
      for (var i = 0; i < names.length; i++) {
        var j = header.indexOf(names[i]);
        if (j >= 0) return j;
      }
      return -1;
    }
    var cPeriod = col(["period", "month", "slot", "date"]);
    var cLabel = col(["label", "y", "malware", "target"]);
    var cScore = col(["score", "logit", "z"]);
    var cProb = col(["prob", "probability", "p"]);
    var cRole = col(["role", "split"]);
    if (cPeriod < 0) throw new Error('Column "period" not found.');
    if (cLabel < 0) throw new Error('Column "label" not found.');
    if (cScore < 0 && cProb < 0) throw new Error('Neither a "score" nor a "prob" column was found.');
    var cValue = cScore >= 0 ? cScore : cProb, valueName = header[cValue];
    var rows = [], invalid = [];
    for (var r = 1; r < lines.length; r++) {
      var f = split(lines[r]);
      var period = f[cPeriod] || "", label = f[cLabel] || "", value = f[cValue] || "";
      var reason = "";
      if (period === "") reason = "empty " + header[cPeriod];
      else if (label === "") reason = "empty " + header[cLabel];
      else if (!NUMBER.test(label) || !(Number(label) === 0 || Number(label) === 1)) {
        reason = header[cLabel] + " is not 0 or 1";
      } else if (value === "") reason = "empty " + valueName;
      else if (!NUMBER.test(value) || !isFinite(Number(value))) reason = valueName + " is not a number";
      if (reason) { invalid.push({ line: numbered[r].line, reason: reason }); continue; }
      rows.push({ period: period, y: Number(label), v: Number(value),
                  role: cRole >= 0 ? (f[cRole] || "").toLowerCase() : "" });
    }
    if (!rows.length) throw new Error("No valid data rows were found.");
    return { rows: rows, skipped: invalid.length, invalid: invalid,
             valueType: cScore >= 0 ? "score" : "prob", hasRole: cRole >= 0 };
  }

  function normalizeRole(s) {
    if (/^c/.test(s)) return "calibration";
    if (/^dev|^val/.test(s)) return "development";
    if (/^dep|^lock|^test|^prod/.test(s)) return "deployment";
    return "";
  }

  /* Group parsed rows into periods in order of first appearance after sorting labels. */
  function groupPeriods(parsed) {
    var map = {}, order = [];
    parsed.rows.forEach(function (r) {
      if (!map[r.period]) { map[r.period] = { label: r.period, v: [], y: [], roles: {} }; order.push(r.period); }
      map[r.period].v.push(r.v);
      map[r.period].y.push(r.y);
      var role = normalizeRole(r.role);
      if (role) map[r.period].roles[role] = (map[r.period].roles[role] || 0) + 1;
    });
    order.sort();
    return order.map(function (k) {
      var p = map[k], roles = Object.keys(p.roles);
      roles.sort(function (a, b) { return p.roles[b] - p.roles[a]; });
      return { label: p.label, v: Float64Array.from(p.v), y: Uint8Array.from(p.y), role: roles[0] || "" };
    });
  }

  /*
   * Periods whose role breaks the order calibration, development, deployment, given the
   * periods in time order. A period conflicts when an earlier period has a later role;
   * periods without one of the three roles are passed over.
   */
  function orderConflicts(periods) {
    var conflicts = [], latest = null;
    periods.forEach(function (p) {
      if (!(p.role in ROLE_ORDER)) return;
      if (latest && ROLE_ORDER[p.role] < ROLE_ORDER[latest.role]) {
        conflicts.push({ label: p.label, role: p.role, after: latest.label, afterRole: latest.role });
      } else {
        latest = p;
      }
    });
    return conflicts;
  }

  /*
   * Build the analysis input from periods with assigned roles, given in time order.
   * Fails when a period is assigned to an earlier role than a period before it, so that
   * no map is fitted or chosen on data later than the periods it is evaluated on.
   */
  function buildData(periods, valueType) {
    var conflicts = orderConflicts(periods);
    if (conflicts.length) {
      throw new Error(
        "Roles must follow the order calibration, development, deployment in time. " +
        conflicts.map(function (c) {
          return "Period " + c.label + " is assigned to " + c.role + " but comes after " +
                 c.afterRole + " period " + c.after + ".";
        }).join(" ") +
        " Periods are ordered by their labels, so labels must sort by time, as 2024-03 does.");
    }
    function toZ(v) {
      if (valueType === "score") return v;
      var z = new Float64Array(v.length);
      for (var i = 0; i < v.length; i++) {
        if (v[i] < 0 || v[i] > 1) throw new Error("Probabilities must lie between 0 and 1.");
        z[i] = logit(clip(v[i]));
      }
      return z;
    }
    var cz = [], cy = [], dev = [], dep = [];
    periods.forEach(function (p) {
      var slot = { label: p.label, z: toZ(p.v), y: p.y };
      if (p.role === "calibration") { cz.push(slot.z); cy.push(slot.y); }
      else if (p.role === "development") dev.push(slot);
      else if (p.role === "deployment") dep.push(slot);
    });
    function concat(parts, Type) {
      var n = parts.reduce(function (s, a) { return s + a.length; }, 0), out = new Type(n), o = 0;
      parts.forEach(function (a) { out.set(a, o); o += a.length; });
      return out;
    }
    if (!cz.length) throw new Error("Assign at least one period to the calibration window.");
    return { calib: { z: concat(cz, Float64Array), y: concat(cy, Uint8Array) }, dev: dev, dep: dep };
  }

  return {
    PCLIP: PCLIP, FAMILIES: FAMILIES, FAMILY_LABEL: FAMILY_LABEL,
    clip: clip, expit: expit, logit: logit, quantile: quantile, binIndex: binIndex,
    makeBinning: makeBinning, dsce: dsce, evidence: evidence, auc: auc, prevalence: prevalence,
    fitLogistic: fitLogistic, fitTemperature: fitTemperature, fitMap: fitMap,
    interceptSlope: interceptSlope, aut: aut, runCadence: runCadence, analyze: analyze,
    decide: decide, parseCsv: parseCsv, groupPeriods: groupPeriods,
    orderConflicts: orderConflicts, buildData: buildData, fmt: fmt
  };
});
