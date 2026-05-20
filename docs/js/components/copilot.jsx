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

// localStorage key for persisted FAB position
var COPILOT_POS_KEY = "cgbw-copilot-pos";

// Read persisted position safely (handles SSR / disabled storage / corrupt JSON)
function loadCopilotPos() {
  try {
    var s = window.localStorage.getItem(COPILOT_POS_KEY);
    if (!s) return null;
    var p = JSON.parse(s);
    if (typeof p.x === "number" && typeof p.y === "number") return p;
  } catch (e) { /* ignore */ }
  return null;
}

// FAB size in CSS px — mirrors the CSS responsive breakpoint
function fabSize() {
  return window.matchMedia && window.matchMedia("(max-width: 720px)").matches ? 46 : 52;
}

// Clamp a (x,y) inside the visible viewport with an 8px margin
function clampPos(p) {
  var s = fabSize();
  var maxX = (window.innerWidth || 0) - s - 8;
  var maxY = (window.innerHeight || 0) - s - 8;
  return {
    x: Math.max(8, Math.min(maxX, p.x)),
    y: Math.max(8, Math.min(maxY, p.y))
  };
}

function Copilot({ onNavigate }) {
  const [doc, setDoc] = React.useState(null);
  const [open, setOpen] = React.useState(false);
  const [input, setInput] = React.useState("");
  const [thread, setThread] = React.useState([]);    // {role:'user'|'bot', text, action?}
  const [typing, setTyping] = React.useState(false);
  // Custom position from drag — null = use default CSS anchor (bottom-left)
  const [pos, setPos] = React.useState(loadCopilotPos);
  const [dragging, setDragging] = React.useState(false);
  const threadEnd = React.useRef(null);
  const inputRef = React.useRef(null);
  const fabRef = React.useRef(null);
  const dragInfo = React.useRef(null);      // { startX, startY, originX, originY, pointerId }
  const wasDragging = React.useRef(false);  // set during pointermove, read in click to suppress

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

  // ─── Persist position whenever it changes (after a drag) ─────────
  React.useEffect(() => {
    if (!pos) return;
    try { window.localStorage.setItem(COPILOT_POS_KEY, JSON.stringify(pos)); }
    catch (e) { /* storage full / disabled — silently ignore */ }
  }, [pos]);

  // ─── Re-clamp position on viewport resize / orientation change ───
  React.useEffect(() => {
    function onResize() {
      setPos(function (p) { return p ? clampPos(p) : p; });
    }
    window.addEventListener("resize", onResize);
    window.addEventListener("orientationchange", onResize);
    return function () {
      window.removeEventListener("resize", onResize);
      window.removeEventListener("orientationchange", onResize);
    };
  }, []);

  // ─── Pointer handlers (drag) ────────────────────────────────────
  // Pointer Events API is supported in Edge 12+, Chrome 55+, Firefox 59+,
  // Safari 13+ — covers every browser the rest of the site supports.
  function onFabPointerDown(e) {
    if (!fabRef.current) return;
    const r = fabRef.current.getBoundingClientRect();
    dragInfo.current = {
      startX: e.clientX, startY: e.clientY,
      originX: r.left,  originY: r.top,
      pointerId: e.pointerId,
    };
    wasDragging.current = false;
    try { fabRef.current.setPointerCapture(e.pointerId); } catch (_) {}
  }

  function onFabPointerMove(e) {
    const info = dragInfo.current;
    if (!info || info.pointerId !== e.pointerId) return;
    const dx = e.clientX - info.startX;
    const dy = e.clientY - info.startY;
    // 5 px dead-zone so a normal click doesn't trip the drag
    if (!wasDragging.current && Math.hypot(dx, dy) < 5) return;
    wasDragging.current = true;
    if (!dragging) setDragging(true);
    setPos(clampPos({ x: info.originX + dx, y: info.originY + dy }));
  }

  function onFabPointerUp(e) {
    const info = dragInfo.current;
    if (!info || info.pointerId !== e.pointerId) return;
    try { fabRef.current.releasePointerCapture(e.pointerId); } catch (_) {}
    dragInfo.current = null;
    if (wasDragging.current) setDragging(false);
    // wasDragging stays true so the synthetic click that fires next is suppressed.
    // It's cleared in onFabClick.
  }

  function onFabClick(e) {
    if (wasDragging.current) {
      wasDragging.current = false;
      e.preventDefault();
      e.stopPropagation();
      return;
    }
    setOpen(function (v) { return !v; });
  }

  // ─── Compute panel style when FAB has been dragged ──────────────
  // Panel opens toward whichever half of the viewport the FAB is in
  // (i.e. into the larger empty area), so it doesn't get pushed off-screen.
  function computePanelStyle() {
    if (!pos) return null;
    const winW = window.innerWidth;
    const winH = window.innerHeight;
    const sz = fabSize();
    const gap = 12;
    const style = {};
    if (pos.x + sz / 2 < winW / 2) {
      style.left  = pos.x + "px";
      style.right = "auto";
    } else {
      style.right = (winW - pos.x - sz) + "px";
      style.left  = "auto";
    }
    if (pos.y + sz / 2 > winH / 2) {
      // FAB in lower half → open panel upward
      style.bottom = (winH - pos.y + gap) + "px";
      style.top    = "auto";
    } else {
      // FAB in upper half → open panel downward
      style.top    = (pos.y + sz + gap) + "px";
      style.bottom = "auto";
    }
    return style;
  }

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
  const fabStyle = pos
    ? { left: pos.x + "px", top: pos.y + "px", right: "auto", bottom: "auto" }
    : null;
  const panelStyle = open ? computePanelStyle() : null;

  return (
    <>
      <button
        ref={fabRef}
        className={"copilot-fab" + (open ? " is-open" : "") + (dragging ? " is-dragging" : "") + (pos ? " is-placed" : "")}
        style={fabStyle}
        onPointerDown={onFabPointerDown}
        onPointerMove={onFabPointerMove}
        onPointerUp={onFabPointerUp}
        onPointerCancel={onFabPointerUp}
        onClick={onFabClick}
        aria-label={open ? "Close Copilot" : "Open Copilot (drag to reposition)"}
        aria-expanded={open}
        title="Click to open · Drag to move">
        {open ? (
          <svg viewBox="0 0 24 24" width="22" height="22" aria-hidden="true">
            <path d="M6 6l12 12M18 6L6 18" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
          </svg>
        ) : (
          <CopilotIcon />
        )}
        {!open && !dragging && <span className="copilot-fab-pulse" aria-hidden="true" />}
      </button>

      <div className={`copilot-panel ${open ? "is-open" : ""}`}
           style={panelStyle}
           role="dialog"
           aria-label="Copilot">
        <header className="copilot-head">
          <span className="copilot-head-mark"><CopilotIcon small /></span>
          <div className="copilot-head-text">
            <span className="copilot-head-name">{doc ? doc.assistant_name : "Copilot"}</span>
            <span className="copilot-head-meta"><span className="copilot-dot" /> demo · offline</span>
          </div>
          {pos && (
            <button className="copilot-head-btn"
                    onClick={() => { setPos(null); try { window.localStorage.removeItem(COPILOT_POS_KEY); } catch(_){} }}
                    title="Move back to default corner"
                    aria-label="Reset position">⤴</button>
          )}
          <button className="copilot-head-btn" onClick={reset} title="New conversation" aria-label="Reset conversation">↺</button>
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

// ─── Brand-matched icon (drop with chat dots, orange Copilot theme) ─
function CopilotIcon({ small }) {
  const sz = small ? 18 : 22;
  return (
    <svg viewBox="0 0 32 32" width={sz} height={sz} aria-hidden="true">
      <defs>
        <linearGradient id="cp-drop" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0%"   stopColor="#f7c98b"/>
          <stop offset="55%"  stopColor="#d8a05a"/>
          <stop offset="100%" stopColor="#8b5e2b"/>
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
