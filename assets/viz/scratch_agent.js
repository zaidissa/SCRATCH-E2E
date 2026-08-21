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

  const ENDPOINTS = CFG.endpoints || [
    { url: "http://localhost:11434/api/chat", kind: "ollama" },
    { url: "http://127.0.0.1:11434/api/chat", kind: "ollama" },
    { url: "http://localhost:8000/v1/chat/completions", kind: "openai" },
  ];

  let llm = null; // { url, kind, model }

  /* ---------------------------------------------------------------------
     Data access
     ------------------------------------------------------------------ */

  const T = (name) => DATA.tables[name] || [];
  const firstTable = () => Object.keys(DATA.tables)[0];
  const numericCols = (rows) =>
    !rows.length ? [] : Object.keys(rows[0]).filter(
      (k) => rows.every((r) => r[k] === null || r[k] === "" || !isNaN(parseFloat(r[k])))
        && rows.some((r) => !isNaN(parseFloat(r[k]))));
  const textCols = (rows) =>
    !rows.length ? [] : Object.keys(rows[0]).filter((k) => !numericCols(rows).includes(k));

  const num = (v) => { const x = parseFloat(v); return isNaN(x) ? null : x; };
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
  add("problems", (q) => /\b(fail|failed|problem|wrong|bad|issue|concern|exclude|drop)\b/.test(q),
    () => {
      const bad = (DATA.findings || []).filter((x) => x.severity === "critical" || x.severity === "warning");
      if (!bad.length) return { text: "No critical or warning findings were registered for this stage.", facts: [] };
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

  // -- distribution of a metric ------------------------------------------
  add("distribution", (q) => /\b(distribution|spread|range|median|mean|average|typical|how many|count)\b/.test(q),
    (q) => {
      const tn = pickTable(q); const rows = T(tn); if (!rows.length) return null;
      const col = pickColumn(q, numericCols(rows));
      if (!col) return { text: `\`${tn}\` has **${rows.length} row(s)**. Columns: ${Object.keys(rows[0]).map((c) => "`" + c + "`").join(", ")}.`, facts: [] };
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

  function pickTable(q) {
    for (const tn of Object.keys(DATA.tables)) if (q.includes(tn.toLowerCase())) return tn;
    return firstTable();
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
    return cols[0] || null;
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
    const tname = pickTable(q);
    const rows = T(tname);
    if (!rows.length) return null;

    const nums = numericCols(rows), txts = textCols(rows);
    const named = Object.keys(rows[0]).filter(
      (c) => q.includes(" " + c.toLowerCase() + " ") || q.includes(c.toLowerCase()));

    const wantHist = /\bhistogram|\bdistribution\b/.test(q);
    const wantScatter = /\bscatter\b|\bvs\b|\bagainst\b/.test(q);

    const numNamed = named.filter((c) => nums.includes(c));
    const txtNamed = named.filter((c) => txts.includes(c));

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
    };
    launch.onclick = () => setOpen(true);

    host.innerHTML = `
      <div class="sxa">
        <div class="sxa-head">
          <span class="sxa-title">${(DATA.meta && DATA.meta.stage) ? DATA.meta.stage + " — ask about this report" : "Ask about this report"}</span>
          <span class="sxa-mode" id="sxa-mode">offline · deterministic</span>
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
      const node = say("bot", md(det.text));
      if (llm && det.intent !== "unknown") {
        const tag = document.createElement("span");
        tag.className = "sxa-tag"; tag.textContent = "phrasing…";
        node.appendChild(tag);
        const better = await phraseWithLLM(q, det);
        if (better) { node.innerHTML = md(better) + '<span class="sxa-tag">local model · grounded in computed facts</span>'; }
        else tag.remove();
      }
    });

    detectLLM().then((f) => {
      llm = f;
      const m = host.querySelector("#sxa-mode");
      if (f) { m.textContent = `offline · local model (${f.model})`; m.classList.add("sxa-live"); }
    });

    // Open on load, after a beat so it slides in over a settled page rather than
    // fighting the widgets for first paint. The summary is the point of the
    // panel; a reader should not have to know to go looking for it. Closing it
    // is one click and the launcher brings it back.
    setTimeout(() => setOpen(true), 400);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", mount);
  else mount();
})();
