// feedback.jsx — Contact page. The address is base64-encoded and only
// revealed when the user clicks the protected mailto chip. No form/box;
// only a single clickable email link + interactive cards.

const TO_B64 = "YWJlZGVtZUBnbWFpbC5jb20=";
const decodeTo = () => { try { return atob(TO_B64); } catch (_) { return ""; } };
const EMAIL_SUBJECT = "CropGBWater – Feedback or Collaboration Request";

function FeedbackPage({ data }) {
  return (
    <section className="feedback">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 07 / Feedback · Collaboration</span>
        </div>
        <h2 className="module-title">
          Reach the creators of the agricultural green- and blue water use Atlas.
        </h2>
        <p className="module-sub">
          We are looking forward to your feedback on the Atlas, scientific
          collaboration, data inquiries, educational use, or partnership
          opportunities.
        </p>
        <ProtectedEmailLink subject={EMAIL_SUBJECT} />
        <WantToHear />
      </div>

      <div className="feedback-links-row">
        <div className="feedback-card">
          <h4>Where else to reach the project</h4>
          <ul className="feedback-links">
            <li><span className="action-key">Paper</span>
                <a href={data?.meta?.paper_url} target="_blank" rel="noopener">{data?.meta?.paper_doi}</a></li>
            <li><span className="action-key">Data</span>
                <a href={data?.meta?.data_url} target="_blank" rel="noopener">Zenodo · {data?.meta?.data_doi}</a></li>
            <li><span className="action-key">Code</span>
                <a href="https://github.com/AbebeDChukalla/CropGBWater" target="_blank" rel="noopener">github.com/AbebeDChukalla/CropGBWater</a></li>
          </ul>
        </div>
      </div>
    </section>
  );
}

// ─── Protected email link ────────────────────────────────────────────
// The address is base64-encoded at rest and never appears as readable text
// on the page until the user actively clicks. First click opens the user's
// email client with the recipient pre-filled; the chip itself only shows
// "Email" as its label. A small lock glyph signals the protection.
function ProtectedEmailLink({ subject }) {
  const [opened, setOpened] = React.useState(false);
  const handleClick = (e) => {
    const to = decodeTo();
    if (!to) return;
    e.preventDefault();
    setOpened(true);
    window.location.href = `mailto:${to}?subject=${encodeURIComponent(subject)}`;
  };
  return (
    <p className="protected-email">
      <span className="protected-email-label">Email:</span>
      <a className="protected-email-link"
         href="#contact"
         rel="noopener"
         aria-label="Open your email client to message the atlas creator"
         onClick={handleClick}>
        <span className="protected-email-lock" aria-hidden="true">🔒</span>
        {opened ? "Opening email client…" : "Click to open in your email client"}
      </a>
    </p>
  );
}

// ─── "What we'd love to hear" — interactive expandable cards ────────
const WTH_ITEMS = [
  { icon: "👍", title: "What works, what doesn't, what's missing",
    body: "Tell us where the Atlas serves your work well — and where it falls short. Concrete examples are gold." },
  { icon: "🗺️", title: "Specific country or crop views you'd like added",
    body: "If there's a country, region, crop or grouping you want to drill into that we don't yet surface, tell us which one and what story you'd want to tell." },
  { icon: "🔍", title: "Discrepancies between the dashboard and your analysis",
    body: "If a number you see here doesn't match your own dataset or paper, we want to know. Send the values + the source." },
  { icon: "📚", title: "Citation or attribution requests",
    body: "Using the Atlas in a publication, course, or report? Happy to help with the right wording, version, and DOI." },
  { icon: "🤝", title: "Collaboration on extending the model",
    body: "New years, new crops, new spatial units, regional deep-dives — open to working together on extensions." },
  { icon: "🎓", title: "Educational use in courses or workshops",
    body: "Using the Atlas to teach water-food-system interactions? Tell us the context and we'll help with material." },
];

function WantToHear() {
  const [open, setOpen] = React.useState(null);
  return (
    <div className="wth">
      <h3 className="wth-title">What we&rsquo;d love to hear</h3>
      <div className="wth-grid">
        {WTH_ITEMS.map((it, i) => {
          const isOpen = open === i;
          return (
            <button key={i}
                    className={`wth-card ${isOpen ? "on" : ""}`}
                    style={{ animationDelay: `${i * 50}ms` }}
                    onClick={() => setOpen(isOpen ? null : i)}
                    aria-expanded={isOpen}>
              <span className="wth-icon" aria-hidden="true">{it.icon}</span>
              <span className="wth-card-body">
                <span className="wth-card-title">{it.title}</span>
                <span className="wth-card-text">{it.body}</span>
              </span>
              <span className="wth-toggle" aria-hidden="true">{isOpen ? "−" : "+"}</span>
            </button>
          );
        })}
      </div>
    </div>
  );
}

window.FeedbackPage = FeedbackPage;
