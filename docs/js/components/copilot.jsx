// copilot.jsx — floating bottom-right chatbot widget.
//
// Honest framing: zero network calls beyond loading data/ai_qa.json. User
// input is matched against `keywords` arrays in the Q&A index; the best
// scoring entry's answer is "typed" out one chunk at a time. If no entry
// scores above 0, the fallback `no_match` text is returned.
//
// Public events:
//   window.dispatchEvent(new CustomEvent("cgbw:open-copilot")) — opens the panel
//
// The widget also exposes a small navigation callback so answers can deep-link.

function Copilot({ onNavigate }) {
  const [doc, setDoc] = React.useState(null);
  const [open, setOpen] = React.useState(false);
  const [input, setInput] = React.useState("");
  const [thread, setThread] = React.useState([]);    // {role:'user'|'bot', text, action?}
  const [typing, setTyping] = React.useState(false);
  const threadEnd = React.useRef(null);
  const inputRef = React.useRef(null);

  // ─── Load Q&A index ──────────────────────────────────────────────
  React.useEffect(() => {
    fetch("data/ai_qa.json")
      .then((r) => r.ok ? r.json() : null)
      .then((d) => setDoc(d))
      .catch(() => setDoc(null));
  }, []);

  // ─── Global open event (other components can open the copilot) ──
  React.useEffect(() => {
    const h = () => setOpen(true);
    window.addEventListener("cgbw:open-copilot", h);
    return () => window.removeEventListener("cgbw:open-copilot", h);
  }, []);

  // ─── Welcome greeting once doc is loaded & panel first opens ────
  React.useEffect(() => {
    if (open && doc && thread.length === 0) {
      setThread([{ role: "bot", text: doc.welcome }]);
      setTimeout(() => inputRef.current && inputRef.current.focus(), 50);
    }
  }, [open, doc, thread.length]);

  // ─── Autoscroll on new messages ─────────────────────────────────
  React.useEffect(() => {
    if (threadEnd.current) threadEnd.current.scrollIntoView({ behavior: "smooth", block: "end" });
  }, [thread, typing]);

  // ─── Send user message ──────────────────────────────────────────
  function send(text) {
    const trimmed = (text || "").trim();
    if (!trimmed || !doc) return;
    setInput("");
    setThread((t) => [...t, { role: "user", text: trimmed }]);
    setTyping(true);

    // Realistic typing delay: 350–800ms depending on length
    const match = pickEntry(trimmed, doc.entries);
    const reply = match ? match.a : doc.no_match;
    const delay = Math.min(900, 350 + Math.round(reply.length / 12));

    setTimeout(() => {
      setThread((t) => [...t, {
        role: "bot",
        text: reply,
        action: match && match.navigates_to ? { target: match.navigates_to } : null
      }]);
      setTyping(false);
    }, delay);
  }

  function handleKey(e) {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      send(input);
    }
  }

  function reset() {
    setThread(doc ? [{ role: "bot", text: doc.welcome }] : []);
    setInput("");
    setTyping(false);
  }

  // ─── Render ─────────────────────────────────────────────────────
  return (
    <>
      <button
        className={`copilot-fab ${open ? "is-open" : ""}`}
        onClick={() => setOpen((v) => !v)}
        aria-label={open ? "Close assistant" : "Open assistant"}
        aria-expanded={open}>
        {open ? (
          <svg viewBox="0 0 24 24" width="22" height="22" aria-hidden="true">
            <path d="M6 6l12 12M18 6L6 18" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
          </svg>
        ) : (
          <CopilotIcon />
        )}
        {!open && <span className="copilot-fab-pulse" aria-hidden="true" />}
      </button>

      <div className={`copilot-panel ${open ? "is-open" : ""}`}
           role="dialog"
           aria-label="CropGBWater Copilot">
        <header className="copilot-head">
          <span className="copilot-head-mark"><CopilotIcon small /></span>
          <div className="copilot-head-text">
            <span className="copilot-head-name">{doc ? doc.assistant_name : "CropGBWater Copilot"}</span>
            <span className="copilot-head-meta"><span className="copilot-dot" /> demo · offline</span>
          </div>
          <button className="copilot-head-btn" onClick={reset} title="New conversation" aria-label="Reset">↺</button>
          <button className="copilot-head-btn" onClick={() => setOpen(false)} aria-label="Close">×</button>
        </header>

        <div className="copilot-body">
          {thread.map((m, i) => (
            <Message key={i} m={m} onNavigate={onNavigate} onClose={() => setOpen(false)} />
          ))}
          {typing && (
            <div className="copilot-msg bot">
              <div className="copilot-bubble typing">
                <span /><span /><span />
              </div>
            </div>
          )}
          <div ref={threadEnd} />
        </div>

        {doc && thread.length <= 1 && (
          <div className="copilot-suggestions" role="group" aria-label="Suggested prompts">
            {doc.suggested_prompts.map((p, i) => (
              <button key={i} className="copilot-chip" onClick={() => send(p)}>{p}</button>
            ))}
          </div>
        )}

        <div className="copilot-input">
          <textarea
            ref={inputRef}
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyDown={handleKey}
            rows={1}
            placeholder={doc ? "Ask about crops, countries, trends…" : "Loading…"}
            aria-label="Ask the assistant" />
          <button className="copilot-send"
                  onClick={() => send(input)}
                  disabled={!input.trim() || typing}
                  aria-label="Send">
            <svg viewBox="0 0 24 24" width="18" height="18" aria-hidden="true">
              <path d="M3 12l18-9-7 18-2-7-9-2z" fill="currentColor" />
            </svg>
          </button>
        </div>
        <p className="copilot-foot">
          Pre-computed demo answers · No live model · <button className="linklike" onClick={() => { setOpen(false); onNavigate && onNavigate("roadmap"); }}>Roadmap</button>
        </p>
      </div>
    </>
  );
}

// ─── Message bubble ────────────────────────────────────────────────
function Message({ m, onNavigate, onClose }) {
  const html = m.role === "bot" ? renderMarkdownish(m.text) : null;
  return (
    <div className={`copilot-msg ${m.role}`}>
      <div className="copilot-bubble">
        {m.role === "bot"
          ? <span dangerouslySetInnerHTML={{ __html: html }} />
          : m.text}
      </div>
      {m.role === "bot" && m.action && (
        <button className="copilot-jump"
                onClick={() => { onClose && onClose(); onNavigate && onNavigate(m.action.target); }}>
          Open {m.action.target} →
        </button>
      )}
    </div>
  );
}

// ─── Tiny markdown subset: **bold**, line breaks, simple "- " bullets, "1." lists ─
function renderMarkdownish(src) {
  const esc = (s) => s.replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  let s = esc(src);
  s = s.replace(/\*\*([^*]+)\*\*/g, "<b>$1</b>");
  // Group consecutive "- " lines into <ul>, "1." lines into <ol>
  const lines = s.split("\n");
  const out = [];
  let listType = null;
  for (const ln of lines) {
    const ulMatch = /^\s*-\s+(.*)$/.exec(ln);
    const olMatch = /^\s*\d+\.\s+(.*)$/.exec(ln);
    if (ulMatch) {
      if (listType !== "ul") { if (listType) out.push(`</${listType}>`); out.push("<ul>"); listType = "ul"; }
      out.push("<li>" + ulMatch[1] + "</li>");
    } else if (olMatch) {
      if (listType !== "ol") { if (listType) out.push(`</${listType}>`); out.push("<ol>"); listType = "ol"; }
      out.push("<li>" + olMatch[1] + "</li>");
    } else {
      if (listType) { out.push(`</${listType}>`); listType = null; }
      if (ln.trim() === "") out.push("<br/>");
      else out.push("<span>" + ln + "</span><br/>");
    }
  }
  if (listType) out.push(`</${listType}>`);
  return out.join("");
}

// ─── Keyword matcher ───────────────────────────────────────────────
// Score = sum over each entry's keywords of: (number of times that keyword
// appears as a substring of the normalised input). Multi-word keywords match
// as a whole phrase. Highest-scoring entry wins; ties broken by first.
function pickEntry(input, entries) {
  if (!entries || entries.length === 0) return null;
  const norm = " " + input.toLowerCase().replace(/[^a-z0-9 ]+/g, " ").replace(/\s+/g, " ") + " ";
  let best = null;
  let bestScore = 0;
  for (const e of entries) {
    let score = 0;
    for (const kw of e.keywords) {
      const k = " " + kw.toLowerCase() + " ";
      // Substring search; multi-word phrases get a length bonus so e.g.
      // "india vs china" outscores entries that only match "india".
      let idx = norm.indexOf(k);
      while (idx !== -1) {
        score += 1 + Math.max(0, kw.split(/\s+/).length - 1);
        idx = norm.indexOf(k, idx + 1);
      }
    }
    if (score > bestScore) { bestScore = score; best = e; }
  }
  return best;
}

// ─── Brand-matched icon (drop with chat dots) ──────────────────────
function CopilotIcon({ small }) {
  const sz = small ? 18 : 22;
  return (
    <svg viewBox="0 0 32 32" width={sz} height={sz} aria-hidden="true">
      <defs>
        <linearGradient id="cp-drop" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0%" stopColor="#7ea36b"/>
          <stop offset="50%" stopColor="#2d5a3d"/>
          <stop offset="100%" stopColor="#1e4a7a"/>
        </linearGradient>
      </defs>
      <path d="M16 3 C 10 11, 8 17, 16 27 C 24 17, 22 11, 16 3 Z"
            fill="url(#cp-drop)" stroke="rgba(0,0,0,0.18)" strokeWidth="0.6"/>
      <circle cx="12" cy="18" r="1.4" fill="#fff"/>
      <circle cx="16" cy="18" r="1.4" fill="#fff"/>
      <circle cx="20" cy="18" r="1.4" fill="#fff"/>
    </svg>
  );
}

window.Copilot = Copilot;
