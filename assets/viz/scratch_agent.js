/* =============================================================================
   scratch_agent.js — offline interpretation agent for SCRATCH-E2E reports
   =============================================================================

   Embedded in the report HTML. Answers questions about the run with NO network:
   the facts, tables and precomputed statistics travel inside the page.

   Two modes, chosen automatically:

     DETERMINISTIC (always available)
       Parses the question, resolves entities (samples, metrics, lineages)
       against the embedded data, and computes the answer — rankings, outliers,
       comparisons, distributions. Every number it states is calculated live
       from the same data the figures were drawn from, so it cannot drift from
       the report. It says what it cannot answer rather than guessing.

     LOCAL MODEL (used when present)
       If a local inference endpoint is reachable — Ollama on 11434, or an
       OpenAI-compatible server such as vLLM/llama.cpp — the agent hands it the
       computed facts as context and lets it phrase the answer. Still no
       internet: the endpoint is on the same host or cluster. If nothing is
       listening, it silently stays deterministic.

   The deterministic path is the product. The model is a phrasing layer on top,
   never the source of a number.
   ========================================================================== */

(function () {
  "use strict";

  const CFG = window.SCRATCH_AGENT_CONFIG || {};
  const DATA = window.SCRATCH_DATA || { findings: [], tables: {}, meta: {} };

  /* ---------------------------------------------------------------------
     Finding schema normalisation
     ---------------------------------------------------------------------
     Two shapes exist. scratch_register() wrote {severity, headline, detail};
     scratch_finding() -- the function every notebook actually calls -- writes
     {type, title, text}. The handlers below were written against the first, so
     a report built from the second answered "6 findings (6 undefined)" and
     printed "undefined — undefined" for every row: the summary was reading keys
     that were not there.
     Normalising once, here, means both shapes work and neither the handlers nor
     the emitting R has to know which produced the page.
     ------------------------------------------------------------------ */
  const SEV_OF_TYPE = { caution: "warning", warning: "warning",
                        important: "critical", danger: "critical",
                        note: "info", tip: "info", info: "info" };
  const SEV_RANK = { critical: 0, warning: 1, good: 2, info: 3 };

  DATA.findings = (DATA.findings || []).map((f) => {
    const sev = f.severity || SEV_OF_TYPE[(f.type || "").toLowerCase()] || "info";
    return {
      ...f,
      severity: sev,
      headline: f.headline || f.title || "(untitled finding)",
      detail:   f.detail   || f.text  || "",
      type:     f.type     || sev,
    };
  }).sort((a, b) => (SEV_RANK[a.severity] ?? 9) - (SEV_RANK[b.severity] ?? 9));

  const ENDPOINTS = CFG.endpoints || [
    { url: "http://localhost:11434/api/chat", kind: "ollama" },
    { url: "http://127.0.0.1:11434/api/chat", kind: "ollama" },
    { url: "http://localhost:8000/v1/chat/completions", kind: "openai" },
  ];

  let llm = null; // { url, kind, model }

  /* ---------------------------------------------------------------------
     Online mode
     ---------------------------------------------------------------------
     The deterministic engine answers from this page. A hosted model adds what
     the page cannot contain: domain knowledge. "Is 11.9% doublets high for 10x?"
     or "what does a kappa of 0.84 mean here?" are questions no embedded table
     can answer.

     THE KEY IS NEVER WRITTEN INTO THE REPORT. These HTML files get emailed,
     copied to USB sticks and posted on shared drives; a key baked into one is a
     key handed to everyone who opens it. It is entered once per browser and kept
     in localStorage, so the report itself stays inert and shareable.

     One thing this cannot do, whatever key you supply: browse the web. The
     reports open over file://, and browsers block cross-origin requests from
     local pages. The model brings its own knowledge; it does not fetch pages.
     ------------------------------------------------------------------ */

  const LS_KEY = "scratch_agent_cloud";

  const CLOUD = {
    anthropic: {
      url: "https://api.anthropic.com/v1/messages",
      model: "claude-sonnet-4-5",
      headers: (k) => ({
        "content-type": "application/json",
        "x-api-key": k,
        "anthropic-version": "2023-06-01",
        // Anthropic blocks browser calls unless this opt-in is present.
        "anthropic-dangerous-direct-browser-access": "true",
      }),
      body: (sys, usr, model) => JSON.stringify({
        model, max_tokens: 900, system: sys,
        messages: [{ role: "user", content: usr }],
      }),
      read: (j) => (j.content && j.content[0] && j.content[0].text) || null,
    },
    openai: {
      url: "https://api.openai.com/v1/chat/completions",
      model: "gpt-4o-mini",
      headers: (k) => ({ "content-type": "application/json", authorization: "Bearer " + k }),
      body: (sys, usr, model) => JSON.stringify({
        model, messages: [{ role: "system", content: sys }, { role: "user", content: usr }],
      }),
      read: (j) => (j.choices && j.choices[0] && j.choices[0].message.content) || null,
    },
  };

  const cloudCfg = () => {
    try { return JSON.parse(localStorage.getItem(LS_KEY) || "null"); } catch (_) { return null; }
  };
  const setCloud = (cfg) => {
    try {
      if (cfg) localStorage.setItem(LS_KEY, JSON.stringify(cfg));
      else localStorage.removeItem(LS_KEY);
    } catch (_) { /* private browsing: online mode is simply unavailable */ }
  };

  // The full report as context. Previously the model was handed only the
  // deterministic engine's answer to reword -- so when that answer was wrong or
  // irrelevant, the model faithfully rephrased something wrong. It now gets the
  // findings and the tables themselves and answers the question directly, with
  // the deterministic result offered as one input among several rather than as
  // the thing to paraphrase.
  function reportContext(maxRows) {
    const t = DATA.tables || {};
    const tables = Object.keys(t).map((name) => {
      const rows = t[name] || [];
      const cols = rows.length ? Object.keys(rows[0]) : [];
      return { name, n_rows: rows.length, columns: cols,
               rows: rows.slice(0, maxRows) };
    });
    return {
      stage: DATA.meta.stage || "unknown",
      project: DATA.meta.project || "unknown",
      findings: (DATA.findings || []).map((f) => ({
        severity: f.severity, title: f.headline, text: strip(f.detail) })),
      tables,
    };
  }

  async function askCloud(question, det) {
    const cfg = cloudCfg();
    if (!cfg || !cfg.key) return null;
    const P = CLOUD[cfg.provider] || CLOUD.anthropic;

    const sys = [
      "You answer questions about a single-cell genomics report, for the scientist who ran it.",
      "",
      "You are given that report's FINDINGS (prose written by the pipeline at render time) and",
      "its TABLES (the actual rows behind the figures). Answer from those.",
      "",
      "Rules:",
      "- Never invent, alter or extrapolate a number. Quote figures exactly as given.",
      "- If the answer is not in the findings or tables, say so plainly and say what IS available.",
      "  Do not substitute a different number that happens to be present.",
      "- Distinguish cells from events, and cells from genes, when the question is ambiguous.",
      "- You MAY add domain knowledge (what a metric means, whether a value is typical, what",
      "  is usually done next). Mark clearly which part is the report and which is your knowledge.",
      "- Be concise. Lead with the answer.",
    ].join("\n");

    const ctx = reportContext(cfg.provider === "anthropic" ? 60 : 30);
    const usr = [
      "QUESTION: " + question,
      "",
      "REPORT CONTEXT (JSON):",
      JSON.stringify(ctx).slice(0, 60000),
      "",
      det && det.intent !== "unknown"
        ? "A deterministic engine also computed this from the same data; use it only if it is relevant:\n" + det.text
        : "(The deterministic engine could not map this question to a computed answer.)",
    ].join("\n");

    try {
      const res = await fetch(P.url, { method: "POST", headers: P.headers(cfg.key),
                                       body: P.body(sys, usr, cfg.model || P.model) });
      if (!res.ok) return { error: `${res.status} ${res.statusText}` };
      const j = await res.json();
      const text = P.read(j);
      return text ? { text } : { error: "empty response" };
    } catch (e) {
      // A browser CORS refusal surfaces here as a TypeError with no detail.
      return { error: e && e.message ? e.message : "request blocked" };
    }
  }

  /* ---------------------------------------------------------------------
     Data access
     ------------------------------------------------------------------ */

  const T = (name) => DATA.tables[name] || [];
  const firstTable = () => Object.keys(DATA.tables)[0];
  // R writes NA into these tables, and a single "NA" used to disqualify an
  // entire column from being numeric -- which is why "plot LLR by state" fell
  // through to an arbitrary column: neutral segments have no LLR, so LLR was
  // classed as text. Missing markers are missing values, not data.
  const MISSING = new Set(["", "na", "nan", "n/a", "null", "none", "inf", "-inf"]);
  const isMissing = (v) => v === null || v === undefined || MISSING.has(String(v).trim().toLowerCase());
  const numericCols = (rows) =>
    !rows.length ? [] : Object.keys(rows[0]).filter(
      (k) => rows.every((r) => isMissing(r[k]) || !isNaN(parseFloat(r[k])))
        && rows.some((r) => !isMissing(r[k]) && !isNaN(parseFloat(r[k]))));
  const textCols = (rows) =>
    !rows.length ? [] : Object.keys(rows[0]).filter((k) => !numericCols(rows).includes(k));

  const num = (v) => { if (isMissing(v)) return null; const x = parseFloat(v); return isNaN(x) ? null : x; };
  const fmt = (x) => (x === null || x === undefined) ? "n/a"
      : (Math.abs(x) >= 1000 ? x.toLocaleString(undefined, { maximumFractionDigits: 0 })
                             : (Math.round(x * 100) / 100).toString());

  function stats(rows, col) {
    const v = rows.map((r) => num(r[col])).filter((x) => x !== null).sort((a, b) => a - b);
    if (!v.length) return null;
    const mean = v.reduce((a, b) => a + b, 0) / v.length;
    const sd = Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / Math.max(v.length - 1, 1));
    const q = (p) => v[Math.min(v.length - 1, Math.floor(p * (v.length - 1)))];
    return { n: v.length, mean, sd, min: v[0], max: v[v.length - 1],
             median: q(0.5), q1: q(0.25), q3: q(0.75) };
  }

  /* ---------------------------------------------------------------------
     Intent handlers — each returns {text, facts} or null
     ------------------------------------------------------------------ */

  const H = [];
  const add = (name, test, run) => H.push({ name, test, run });

  // -- what happened / summarise -----------------------------------------
  add("summary", (q) => /\bsummar\w*|\boverview\b|what happened|explain the run|\bbrief\b|how did .* go/.test(q),
    () => {
      const f = DATA.findings || [];
      const bySev = {};
      f.forEach((x) => (bySev[x.severity] = (bySev[x.severity] || 0) + 1));
      const lines = f.slice(0, 8).map((x) => `- **${x.severity}** — ${x.headline}${x.detail ? ": " + strip(x.detail) : ""}`);
      return {
        text: [`This run registered **${f.length} finding(s)**` +
               (Object.keys(bySev).length ? ` (${Object.entries(bySev).map(([k, v]) => `${v} ${k}`).join(", ")}).` : "."),
               "", ...lines].join("\n"),
        facts: f,
      };
    });

  // -- failures / problems ------------------------------------------------
  // "needs attention" and "worth checking" are how people actually phrase this —
  // and two of them are the panel's own suggestion chips, which fell through to
  // the "I could not map that" fallback until they were added here. A chip that
  // does not work is worse than no chip.
  add("problems", (q) => /\b(fail|failed|problem|wrong|bad|issue|concern|exclude|drop|caution|warning|attention|worry|risk|check)\b/.test(q),
    () => {
      const bad = (DATA.findings || []).filter((x) => x.severity === "critical" || x.severity === "warning");
      if (!bad.length) {
        const n = (DATA.findings || []).length;
        return { text: n
          ? `Nothing is flagged as a caution here. All ${n} finding(s) on this page are informational — ask me to **summarise this report** to see them.`
          : "No findings were registered for this stage.", facts: [] };
      }
      return {
        text: bad.map((x) => `- **${x.severity}** — ${x.headline}${x.detail ? ": " + strip(x.detail) : ""}`).join("\n"),
        facts: bad,
      };
    });

  // -- rare / outlier / unusual ------------------------------------------
  add("outliers", (q) => /\b(rare|outlier|unusual|anomal|odd|strange|stand ?out|atypical|surpris)\b/.test(q),
    (q) => {
      const tn = pickTable(q); const rows = T(tn);
      if (!rows.length) return null;
      const idc = textCols(rows)[0];
      const out = [];
      numericCols(rows).forEach((c) => {
        const s = stats(rows, c); if (!s || s.sd === 0) return;
        rows.forEach((r) => {
          const z = (num(r[c]) - s.mean) / s.sd;
          if (Math.abs(z) >= 2) out.push({ id: r[idc], metric: c, value: num(r[c]), z, median: s.median });
        });
      });
      out.sort((a, b) => Math.abs(b.z) - Math.abs(a.z));
      if (!out.length) return { text: `No value in **${tn}** sits more than 2 SD from its column mean — the cohort is homogeneous on the measured metrics.`, facts: [] };
      return {
        text: [`**${out.length} outlying value(s)** in \`${tn}\` (|z| ≥ 2), strongest first:`, "",
          ...out.slice(0, 8).map((o) =>
            `- **${o.id}** — \`${o.metric}\` = ${fmt(o.value)} (${o.z > 0 ? "+" : ""}${o.z.toFixed(1)} SD; cohort median ${fmt(o.median)})`)].join("\n"),
        facts: out.slice(0, 20),
      };
    });

  // -- top / rank ---------------------------------------------------------
  add("ranking", (q) => /\b(top|highest|lowest|best|worst|rank|most|least|largest|smallest)\b/.test(q),
    (q) => {
      const tn = pickTable(q); const rows = T(tn); if (!rows.length) return null;
      const col = pickColumn(q, numericCols(rows)); if (!col) return null;
      const asc = /\b(lowest|worst|least|smallest|fewest|bottom)\b/.test(q);
      const n = (q.match(/\btop\s+(\d+)|\b(\d+)\s+(highest|lowest|most|least)/) || [])
        .filter(Boolean).map(Number).find((x) => !isNaN(x)) || 5;
      const idc = textCols(rows)[0];
      const sorted = rows.filter((r) => num(r[col]) !== null)
        .sort((a, b) => asc ? num(a[col]) - num(b[col]) : num(b[col]) - num(a[col]));
      const s = stats(rows, col);
      return {
        text: [`**${asc ? "Lowest" : "Highest"} ${Math.min(n, sorted.length)} by \`${col}\`** (cohort median ${fmt(s.median)}):`, "",
          ...sorted.slice(0, n).map((r, i) => `${i + 1}. **${r[idc]}** — ${fmt(num(r[col]))}`)].join("\n"),
        facts: sorted.slice(0, n),
      };
    });

  // -- compare two entities ----------------------------------------------
  add("compare", (q) => /\b(compare|versus|vs\.?|difference between|against)\b/.test(q),
    (q) => {
      const tn = pickTable(q); const rows = T(tn); if (!rows.length) return null;
      const idc = textCols(rows)[0];
      const hits = rows.filter((r) => q.includes(String(r[idc]).toLowerCase()));
      if (hits.length < 2) return null;
      const cols = numericCols(rows);
      const [a, b] = hits;
      return {
        text: [`**${a[idc]}** vs **${b[idc]}**:`, "",
          ...cols.map((c) => {
            const va = num(a[c]), vb = num(b[c]);
            if (va === null || vb === null) return null;
            const d = vb - va, pc = va !== 0 ? (100 * d / Math.abs(va)) : NaN;
            return `- \`${c}\` — ${fmt(va)} → ${fmt(vb)} (${d > 0 ? "+" : ""}${fmt(d)}${isNaN(pc) ? "" : `, ${d > 0 ? "+" : ""}${pc.toFixed(0)}%`})`;
          }).filter(Boolean)].join("\n"),
        facts: hits,
      };
    });

  // -- a specific entity --------------------------------------------------
  // Registered late in the file but PRIORITISED at match time (see the reorder
  // below): naming a sample is a strong signal and must beat generic intents.
  add("entity", (q) => !!findEntity(q),
    (q) => {
      for (const tn of Object.keys(DATA.tables)) {
        const rows = T(tn); if (!rows.length) continue;
        const idc = textCols(rows)[0];
        const hit = rows.find((r) => q.includes(String(r[idc]).toLowerCase()));
        if (hit) {
          return {
            text: [`**${hit[idc]}** — every recorded value:`, "",
              ...Object.entries(hit).filter(([k]) => k !== idc)
                .map(([k, v]) => `- \`${k}\` = ${v === "" ? "n/a" : v}`)].join("\n"),
            facts: [hit],
          };
        }
      }
      return null;
    });

  // -- answer from the findings ------------------------------------------
  // Ranked ahead of the table handlers. Most substantive numbers in these
  // reports are stated in findings ("Numbat called 9,227"), not in a registered
  // table -- so a question naming a caller, a stage or a metric is far more
  // likely to be answerable from prose than from a column. Quoting the finding
  // is also honest: the reader sees the sentence the report itself wrote.
  // Substring matching ranked the wrong finding first: "cnv" matches INSIDE
  // "inferCNV", so a question about Numbat's CNVs scored the inferCNV coverage
  // finding just as highly as the segment count that actually answers it.
  // Whole-word matching separates them.
  const hasTerm = (hay, t) => new RegExp("\\b" + t.replace(/[.*+?^${}()|[\]\\]/g, "\\$&") + "\\b").test(hay);

  add("findings-search", (q) => {
    const f = DATA.findings || [];
    if (!f.length) return false;
    return terms(q).some((t) => f.some((x) => hasTerm((x.headline + " " + x.detail).toLowerCase(), t)));
  }, (q) => {
    const f = DATA.findings || [];
    const ts = terms(q);
    const scored = f.map((x) => {
      const hay = (x.headline + " " + x.detail).toLowerCase();
      // Title matches count double: a finding whose HEADING is about the thing
      // asked about is more likely the answer than one that mentions it in passing.
      const body = ts.filter((t) => hasTerm(hay, t)).length;
      const head = ts.filter((t) => hasTerm(x.headline.toLowerCase(), t)).length;
      return { x, hits: body + head };
    }).filter((r) => r.hits > 0).sort((a, b) => b.hits - a.hits);
    if (!scored.length) return null;

    // Relevance floor. Every finding in a report shares vocabulary -- numbat,
    // cnv, cells, sample -- so "any term matches any finding" fired on nearly
    // every question and returned the SAME two or three findings regardless of
    // what was asked. Require either a title hit or two distinct matching terms,
    // and fall through to the other handlers otherwise.
    const strong = scored.filter((r) =>
      r.hits >= 2 || ts.some((t) => hasTerm(r.x.headline.toLowerCase(), t)));
    if (!strong.length) return null;
    const lines = strong.slice(0, 3).map((r) =>
      `- **${strip(r.x.headline)}** — ${strip(r.x.detail)}`);
    return {
      text: ["From this report's findings:", "", ...lines, "",
             "_If you need a number that is not stated above, it is not in this page's data._"].join("\n"),
      facts: strong.slice(0, 3).map((r) => r.x),
    };
  });

  // -- distribution of a metric ------------------------------------------
  add("distribution", (q) => /\b(distribution|spread|range|median|mean|average|typical|how many|count)\b/.test(q),
    (q) => {
      const tn = pickTable(q); const rows = T(tn); if (!rows.length) return null;
      const col = pickColumn(q, numericCols(rows));
      // Fall through rather than describe an arbitrary column: "how many X" with
      // no recognised column is a question about something this table may not
      // contain, and answering it with whatever was registered first invents an
      // answer.
      if (!col) return null;
      const s = stats(rows, col);
      return {
        text: [`\`${col}\` across ${s.n} entries in \`${tn}\`:`, "",
          `- median **${fmt(s.median)}** (IQR ${fmt(s.q1)}–${fmt(s.q3)})`,
          `- mean ${fmt(s.mean)} ± ${fmt(s.sd)} SD`,
          `- range ${fmt(s.min)} to ${fmt(s.max)}`].join("\n"),
        facts: [s],
      };
    });

  function strip(s) { return String(s).replace(/\*\*/g, ""); }

  // Content words worth matching on. The stop list keeps "how many" and friends
  // from matching every finding in the page.
  const STOP = new Set(["the","a","an","of","in","on","for","to","and","or","is","are",
    "was","were","how","many","much","what","which","this","that","these","those","by",
    "with","from","do","does","did","can","you","me","my","it","its","have","has","any",
    "report","stage","run","cells","cell","show","tell","give"]);
  function terms(q) {
    return String(q).toLowerCase().split(/[^a-z0-9_]+/)
      .filter((w) => w.length > 2 && !STOP.has(w));
  }

  // Does the question name a row in any embedded table?
  function findEntity(q) {
    for (const tn of Object.keys(DATA.tables)) {
      const rows = T(tn); if (!rows.length) continue;
      const idc = textCols(rows)[0];
      const hit = rows.find((r) => String(r[idc]).length > 2 && q.includes(String(r[idc]).toLowerCase()));
      if (hit) return { table: tn, row: hit, idc };
    }
    return null;
  }

  // Intent priority. A named entity or an explicit comparison is far more
  // specific than "summarise", so those are consulted first regardless of the
  // order the handlers happen to be defined in.
  const PRIORITY = ["compare", "ranking", "outliers", "entity", "problems", "distribution", "summary"];
  function orderedHandlers() {
    const byName = Object.fromEntries(H.map((h) => [h.name, h]));
    return PRIORITY.map((n) => byName[n]).filter(Boolean)
      .concat(H.filter((h) => !PRIORITY.includes(h.name)));
  }

  // Returns null when the question names no table. It used to fall back to
  // firstTable(), which is how "how many CNV identified by numbat" ended up
  // being answered with statistics about an unrelated calibration column: a
  // confident, plausible, completely wrong number. Guessing is worse than
  // refusing here, because the reader cannot tell the difference.
  function pickTable(q, { fallback = false } = {}) {
    for (const tn of Object.keys(DATA.tables)) if (q.includes(tn.toLowerCase())) return tn;
    if (!fallback) return null;
    // Explicit fallback (used only for "plot the main table"): the biggest table
    // is the most informative default, not whichever happened to be registered first.
    let best = null, n = -1;
    for (const tn of Object.keys(DATA.tables)) {
      const len = (DATA.tables[tn] || []).length;
      if (len > n) { n = len; best = tn; }
    }
    return best;
  }
  function pickColumn(q, cols) {
    let best = null, bestLen = 0;
    cols.forEach((c) => {
      const c2 = c.toLowerCase().replace(/[_.]/g, " ");
      if ((q.includes(c.toLowerCase()) || q.includes(c2)) && c.length > bestLen) { best = c; bestLen = c.length; }
    });
    if (best) return best;
    const syn = { cell: /\bcell|barcode/, gene: /\bgene|complexity|feature/,
                  mito: /\bmito|mitochond/, doublet: /\bdoublet/,
                  umi: /\bumi|count|depth/ };
    for (const c of cols) for (const [k, re] of Object.entries(syn))
      if (re.test(q) && c.toLowerCase().includes(k)) return c;
    // No match: say so. Returning cols[0] made every unrecognised question look
    // answered.
    return null;
  }

  /* ---------------------------------------------------------------------
     Answer
     ------------------------------------------------------------------ */

  function answerDeterministic(question) {
    const q = " " + question.toLowerCase().trim() + " ";
    for (const h of orderedHandlers()) {
      if (!h.test(q)) continue;
      let r = null;
      try { r = h.run(q); } catch (e) { r = null; }
      if (r && r.text) return { ...r, intent: h.name };
    }
    return {
      intent: "unknown",
      text: ["I could not map that to the data in this report. I can answer:", "",
        "- **what happened** — the run's findings, ranked by severity",
        "- **what failed / what's concerning**",
        "- **rare or unusual patterns** — values ≥ 2 SD from the cohort",
        "- **top / lowest N by a metric**",
        "- **compare A vs B**",
        "- **tell me about \\<sample\\>**",
        "- **distribution of a metric**", "",
        `Available tables: ${Object.keys(DATA.tables).map((t) => "`" + t + "`").join(", ") || "none"}.`].join("\n"),
      facts: [],
    };
  }

  async function detectLLM() {
    for (const e of ENDPOINTS) {
      try {
        const base = e.url.replace(/\/(api\/chat|v1\/chat\/completions)$/, "");
        const probe = e.kind === "ollama" ? base + "/api/tags" : base + "/v1/models";
        const ctl = new AbortController();
        const t = setTimeout(() => ctl.abort(), 700);
        const res = await fetch(probe, { signal: ctl.signal });
        clearTimeout(t);
        if (!res.ok) continue;
        const j = await res.json();
        const model = CFG.model
          || (j.models && j.models[0] && (j.models[0].name || j.models[0].model))
          || (j.data && j.data[0] && j.data[0].id);
        if (model) return { ...e, model };
      } catch (_) { /* nothing listening; stay deterministic */ }
    }
    return null;
  }

  async function phraseWithLLM(question, det) {
    if (!llm) return null;
    const sys = "You explain single-cell genomics QC results. You are given FACTS computed "
      + "from the data. Use ONLY those facts. Never invent a number. Be concise and concrete. "
      + "If the facts do not answer the question, say so plainly.";
    const usr = `Question: ${question}\n\nFACTS (computed from the report data):\n${det.text}\n\n`
      + `Structured facts:\n${JSON.stringify(det.facts).slice(0, 4000)}\n\n`
      + `Stage: ${DATA.meta.stage || "unknown"}. Project: ${DATA.meta.project || "unknown"}.`;
    try {
      const body = llm.kind === "ollama"
        ? { model: llm.model, stream: false, messages: [{ role: "system", content: sys }, { role: "user", content: usr }] }
        : { model: llm.model, messages: [{ role: "system", content: sys }, { role: "user", content: usr }], temperature: 0.2 };
      const ctl = new AbortController();
      const t = setTimeout(() => ctl.abort(), 30000);
      const res = await fetch(llm.url, { method: "POST", headers: { "Content-Type": "application/json" },
                                         body: JSON.stringify(body), signal: ctl.signal });
      clearTimeout(t);
      if (!res.ok) return null;
      const j = await res.json();
      return (j.message && j.message.content) || (j.choices && j.choices[0].message.content) || null;
    } catch (_) { return null; }
  }

  /* ---------------------------------------------------------------------
     UI
     ------------------------------------------------------------------ */

  function md(s) {
    return String(s)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
      .replace(/\*\*(.+?)\*\*/g, "<strong>$1</strong>")
      .replace(/`([^`]+)`/g, "<code>$1</code>")
      .replace(/^- (.*)$/gm, "<li>$1</li>")
      .replace(/(<li>[\s\S]*<\/li>)/, "<ul>$1</ul>")
      .replace(/\n{2,}/g, "<br><br>").replace(/\n/g, "<br>");
  }

  /* ---------------------------------------------------------------------
     Charts the agent draws on request
     ---------------------------------------------------------------------
     Plotly is already in the page (the report's own figures are widgets), so a
     chart costs no new dependency. If a report happens to carry no widget,
     Plotly is absent and these handlers say so instead of throwing.

     Deliberately limited to the registered tables. The agent can plot what the
     report computed; it cannot invent a variable that was never measured.
     ------------------------------------------------------------------ */

  const havePlotly = () => typeof window.Plotly !== "undefined";

  function plotTheme(title) {
    const dark = window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
    return {
      title: { text: title, font: { size: 13 } },
      margin: { l: 52, r: 16, t: 38, b: 46 },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: dark ? "#e6edf3" : "#24292f", size: 11 },
      showlegend: false,
    };
  }

  function drawPlot(node, traces, title) {
    if (!havePlotly()) {
      node.insertAdjacentHTML("beforeend",
        "<br><em>This report carries no Plotly runtime, so I cannot draw here. " +
        "The numbers above are still computed from the data.</em>");
      return;
    }
    const div = document.createElement("div");
    div.className = "sxa-plot";
    node.appendChild(div);
    window.Plotly.newPlot(div, traces, plotTheme(title),
                          { displayModeBar: false, responsive: true });
  }

  // "plot|chart|histogram|scatter ..." -> { traces, title, summary }
  function buildChart(q) {
    // fallback:true here, unlike the text handlers: "plot the main table" is an
    // explicit request to choose one, so choosing the LARGEST is reasonable --
    // and the summary line names it, so the reader can see what was picked.
    const tname = pickTable(q, { fallback: true });
    const rows = T(tname);
    if (!rows.length) return null;

    const nums = numericCols(rows), txts = textCols(rows);
    const named = Object.keys(rows[0]).filter(
      (c) => q.includes(" " + c.toLowerCase() + " ") || q.includes(c.toLowerCase()));

    const wantHist = /\bhistogram|\bdistribution\b/.test(q);
    const wantScatter = /\bscatter\b|\bvs\b|\bagainst\b/.test(q);

    const numNamed = named.filter((c) => nums.includes(c));
    const txtNamed = named.filter((c) => txts.includes(c));

    // "plot X by Y" names the measure first and the grouping second. Reordering
    // the candidate lists to express that turned out to interact badly with the
    // fallback below ("plot llr by state" drew CHROM), so an explicit A-by-B
    // request now builds its spec directly and returns. Less clever, and it does
    // what it says.
    const byM = q.match(/\b([a-z0-9_.]+)\s+(?:by|per|across|against|vs\.?)\s+([a-z0-9_.]+)/);
    if (byM) {
      const find = (w) => Object.keys(rows[0]).find((c) =>
        c.toLowerCase() === w || c.toLowerCase().replace(/[_.]/g, "") === w.replace(/[_.]/g, ""));
      const yc = find(byM[1]), xc = find(byM[2]);
      if (yc && xc && nums.includes(yc)) {
        const bothNum = nums.includes(xc);
        if (bothNum && wantScatter) {
          return {
            title: `${yc} vs ${xc} — ${tname}`,
            summary: `Scatter of **${yc}** against **${xc}** over ${rows.length} rows of \`${tname}\`.`,
            traces: [{ x: rows.map((r) => num(r[xc])), y: rows.map((r) => num(r[yc])),
                       mode: "markers", type: "scattergl",
                       marker: { size: 7, color: "#2a78d6", opacity: .8 } }],
          };
        }
        // Otherwise treat x as a grouping and SUM the measure within it, which is
        // what "size by chromosome" means -- not one bar per row.
        const agg = {};
        rows.forEach((r) => {
          const k = String(r[xc]); const v = num(r[yc]);
          if (v === null) return;
          agg[k] = (agg[k] || 0) + v;
        });
        const keys = Object.keys(agg).sort((a, b) => (parseFloat(a) - parseFloat(b)) || a.localeCompare(b));
        return {
          title: `${yc} by ${xc} — ${tname}`,
          summary: `**${yc}** summed within each **${xc}** (${keys.length} group(s), ${rows.length} rows of \`${tname}\`).`,
          traces: [{ x: keys, y: keys.map((k) => agg[k]), type: "bar",
                     marker: { color: "#2a78d6" } }],
        };
      }
    }

    if (wantScatter && numNamed.length >= 2) {
      const [x, y] = numNamed;
      return {
        title: `${y} vs ${x} — ${tname}`,
        summary: `Scatter of **${y}** against **${x}** over ${rows.length} rows of \`${tname}\`.`,
        traces: [{ x: rows.map((r) => num(r[x])), y: rows.map((r) => num(r[y])),
                   text: rows.map((r) => txts.length ? r[txts[0]] : ""),
                   mode: "markers", type: "scattergl",
                   marker: { size: 7, color: "#2a78d6", opacity: .8 } }],
      };
    }

    if (wantHist && numNamed.length >= 1) {
      const x = numNamed[0];
      return {
        title: `Distribution of ${x} — ${tname}`,
        summary: `Histogram of **${x}** over ${rows.length} rows of \`${tname}\`.`,
        traces: [{ x: rows.map((r) => num(r[x])), type: "histogram",
                   marker: { color: "#2a78d6" } }],
      };
    }

    // Default: a bar of one numeric by one categorical — the shape most of these
    // tables actually are (per-sample, per-cluster, per-cell-type summaries).
    const yv = numNamed[0] || nums[0];
    const xv = txtNamed[0] || txts[0];
    if (!yv || !xv) return null;
    const idx = rows.map((r, i) => i).sort((a, b) => (num(rows[b][yv]) || 0) - (num(rows[a][yv]) || 0)).slice(0, 30);
    return {
      title: `${yv} by ${xv} — ${tname}`,
      summary: `**${yv}** by **${xv}** from \`${tname}\`${idx.length < rows.length ? `, top ${idx.length}` : ""}.`,
      traces: [{ x: idx.map((i) => rows[i][xv]), y: idx.map((i) => num(rows[i][yv])),
                 type: "bar", marker: { color: "#2a78d6" } }],
    };
  }

  function mount() {
    const host = document.getElementById("scratch-agent");
    if (!host) return;

    // Launcher lives outside the panel so it stays reachable when the panel is
    // closed. The panel opens on first load: a chat box nobody notices is the
    // same as no chat box, which is what the previous inline version was.
    const launch = document.createElement("button");
    launch.className = "sxa-launch";
    launch.innerHTML = '<span class="sxa-dot"></span> Ask about this report';
    document.body.appendChild(launch);

    const setOpen = (open) => {
      host.classList.toggle("sxa-open", open);
      launch.style.display = open ? "none" : "flex";
      // Shift the document instead of covering it. Without this the panel sits
      // on top of Quarto's right-hand TOC.
      document.body.classList.toggle("sxa-shifted", open);
    };
    launch.onclick = () => setOpen(true);

    host.innerHTML = `
      <div class="sxa">
        <div class="sxa-head">
          <span class="sxa-title">${(DATA.meta && DATA.meta.stage) ? DATA.meta.stage + " — ask about this report" : "Ask about this report"}</span>
          <span class="sxa-mode" id="sxa-mode">offline · deterministic</span>
          <button class="sxa-close" id="sxa-cfg" title="Online mode" aria-label="Online mode">&#9881;</button>
          <button class="sxa-close" id="sxa-close" title="Close" aria-label="Close">&times;</button>
        </div>
        <div class="sxa-log" id="sxa-log"></div>
        <div class="sxa-chips" id="sxa-chips"></div>
        <form class="sxa-form" id="sxa-form" autocomplete="off">
          <input id="sxa-in" placeholder="e.g. do you see any rare patterns?" aria-label="Ask about this run">
          <button type="submit">Ask</button>
        </form>
      </div>`;

    const log = host.querySelector("#sxa-log");
    const chips = host.querySelector("#sxa-chips");
    ["Summarise this report", "What needs attention?", "Any rare patterns?",
     "Plot the main table", "Top 5 by cell count"]
      .forEach((s) => {
        const b = document.createElement("button");
        b.className = "sxa-chip"; b.type = "button"; b.textContent = s;
        b.onclick = () => { host.querySelector("#sxa-in").value = s; host.querySelector("#sxa-form").requestSubmit(); };
        chips.appendChild(b);
      });

    function say(who, html, cls) {
      const d = document.createElement("div");
      d.className = "sxa-msg sxa-" + who + (cls ? " " + cls : "");
      d.innerHTML = html;
      log.appendChild(d); log.scrollTop = log.scrollHeight;
      return d;
    }

    host.querySelector("#sxa-close").onclick = () => setOpen(false);

    // ---- online mode -----------------------------------------------------
    const modeEl = host.querySelector("#sxa-mode");
    function refreshMode() {
      const c = cloudCfg();
      if (c && c.key) {
        modeEl.textContent = "online · " + (c.provider || "anthropic");
        modeEl.classList.add("sxa-live");
      } else if (llm) {
        modeEl.textContent = "offline · local model (" + llm.model + ")";
        modeEl.classList.add("sxa-live");
      } else {
        modeEl.textContent = "offline · deterministic";
        modeEl.classList.remove("sxa-live");
      }
    }

    host.querySelector("#sxa-cfg").onclick = () => {
      const cur = cloudCfg() || {};
      if (cur.key) {
        if (confirm("Online mode is on (" + (cur.provider || "anthropic") + ").\n\nTurn it off and forget the key stored in this browser?")) {
          setCloud(null); refreshMode();
          say("bot", md("Online mode off. The key has been removed from this browser. "
            + "Answers come from the deterministic engine and this page's data."));
        }
        return;
      }
      const provider = (prompt("Online mode.\n\nProvider — type 'anthropic' or 'openai':", "anthropic") || "").trim().toLowerCase();
      if (!provider) return;
      if (!CLOUD[provider]) { say("bot", md("Unknown provider. Use `anthropic` or `openai`.")); return; }
      const key = (prompt("Paste your " + provider + " API key.\n\nIt is stored in THIS BROWSER only (localStorage) and is never written into the report file, so sharing this HTML does not share your key.") || "").trim();
      if (!key) return;
      setCloud({ provider, key, model: CLOUD[provider].model });
      refreshMode();
      say("bot", md("**Online mode on** (" + provider + ", " + CLOUD[provider].model + "). "
        + "I can now bring outside knowledge — what a metric means, whether a value is typical, "
        + "what is usually done next — on top of this report's computed numbers. "
        + "I still cannot browse the web: browsers block that from a local file.\n\n"
        + "The key lives in this browser only. Click the gear again to remove it."));
    };

    // Open with the summary already written. The first thing a reader wants is
    // "what does this report say", and making them ask for it is friction.
    const openingSummary = () => {
      const f = DATA.findings || [];
      const sev = (t) => f.filter((x) => (x.type || "").toLowerCase() === t).length;
      const nT = Object.keys(DATA.tables || {}).length;
      const bits = [];
      bits.push(`**${DATA.meta && DATA.meta.stage ? DATA.meta.stage : "This stage"}**`
        + `${DATA.meta && DATA.meta.project ? " · " + DATA.meta.project : ""}`
        + ` — ${f.length} finding(s), ${nT} table(s) in this page.`);
      const caution = f.filter((x) => /caution|important/i.test(x.type || ""));
      if (caution.length) {
        bits.push("", `**${caution.length} need attention:**`);
        caution.slice(0, 4).forEach((x) => bits.push(`- ${strip(x.title || "")}: ${strip(x.text || "").slice(0, 190)}`));
      } else if (f.length) {
        bits.push("", "Nothing flagged as a caution. Headline findings:");
        f.slice(0, 3).forEach((x) => bits.push(`- ${strip(x.title || "")}: ${strip(x.text || "").slice(0, 170)}`));
      }
      bits.push("", "Ask me anything about these numbers, or say _plot \`<column>\` by \`<column>\`_ and I will draw it.");
      return bits.join("\n");
    };
    say("bot", md(openingSummary()));

    host.querySelector("#sxa-form").addEventListener("submit", async (e) => {
      e.preventDefault();
      const inp = host.querySelector("#sxa-in");
      const q = inp.value.trim(); if (!q) return;
      inp.value = "";
      say("user", md(q));

      // Plot requests are handled before the deterministic text handlers: they
      // produce a figure plus a sentence, not a paragraph.
      if (/\bplot\b|\bchart\b|\bgraph\b|\bhistogram\b|\bdraw\b|\bvisuali[sz]e\b/.test(" " + q.toLowerCase() + " ")) {
        const spec = buildChart(" " + q.toLowerCase() + " ");
        if (spec) { drawPlot(say("bot", md(spec.summary)), spec.traces, spec.title); return; }
        say("bot", md("I could not work out what to plot from that. Name a table column — "
          + "e.g. _plot n_observed_cells by sample_id_, or _histogram of percent_mito_. "
          + `Columns I have: ${Object.keys(DATA.tables).map((t) => "`" + t + "`").join(", ") || "none"}.`));
        return;
      }

      const det = answerDeterministic(q);
      const online = cloudCfg() && cloudCfg().key;

      // With online mode on, the model answers the question -- it is not shown a
      // deterministic answer to reword. Showing that answer first was actively
      // misleading when it was the wrong answer, which for open-ended questions
      // it usually is.
      const node = say("bot", md(online ? "_Reading the report…_" : det.text));

      if (online) {
        const tag = document.createElement("span");
        tag.className = "sxa-tag"; tag.textContent = "asking the model…";
        node.appendChild(tag);
        const r = await askCloud(q, det);
        if (r && r.text) {
          node.innerHTML = md(r.text) + '<span class="sxa-tag">online model · numbers from this report</span>';
        } else {
          // Fall back to the offline answer rather than leaving a dead message.
          node.innerHTML = md(det.text)
            + '<span class="sxa-tag">online call failed ('
            + ((r && r.error) || "no response") + ") — computed answer shown</span>";
        }
        return;
      }

      if (llm && det.intent !== "unknown") {
        const tag = document.createElement("span");
        tag.className = "sxa-tag"; tag.textContent = "phrasing…";
        node.appendChild(tag);
        const better = await phraseWithLLM(q, det);
        if (better) { node.innerHTML = md(better) + '<span class="sxa-tag">local model · grounded in computed facts</span>'; }
        else tag.remove();
      }
    });

    detectLLM().then((f) => { llm = f; refreshMode(); });
    refreshMode();

    // Open on load, after a beat so it slides in over a settled page rather than
    // fighting the widgets for first paint. The summary is the point of the
    // panel; a reader should not have to know to go looking for it. Closing it
    // is one click and the launcher brings it back.
    setTimeout(() => setOpen(true), 400);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", mount);
  else mount();
})();
