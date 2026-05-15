// feedback.jsx — Contact page. Email is hidden everywhere on the page;
// the destination is base64-encoded so it doesn't appear in plain text.
// Posts to FormSubmit; on failure (network/CORS), falls back to a native
// form submit (full page POST) which works even when fetch() is blocked.

// Base64 of "abedeme@gmail.com" — decoded only at submit time so the address
// is never present in plain HTML or visible text on the page.
const TO_B64 = "YWJlZGVtZUBnbWFpbC5jb20=";
const decodeTo = () => {
  try { return atob(TO_B64); } catch (_) { return ""; }
};

const EMAIL_SUBJECT = "CropGBWater – Feedback or Collaboration Request";

function FeedbackPage({ data }) {
  const [status, setStatus] = React.useState("idle"); // idle | sending | success | error
  const [error,  setError]  = React.useState(null);
  const formRef = React.useRef(null);

  const fallbackSubmit = () => {
    // Force a regular form POST as a fallback. The form's action is set on
    // demand (so the address is never in the static HTML).
    const f = formRef.current;
    if (!f) return;
    f.action = `https://formsubmit.co/${decodeTo()}`;
    f.method = "POST";
    // Remove the AJAX handler so the browser does a regular submit
    f.onsubmit = null;
    f.submit();
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    const f = e.currentTarget;
    setStatus("sending");
    setError(null);
    try {
      const fd = new FormData(f);
      const res = await fetch(`https://formsubmit.co/ajax/${decodeTo()}`, {
        method: "POST",
        headers: { "Accept": "application/json" },
        body: fd,
      });
      const json = await res.json().catch(() => ({}));
      if (res.ok && (json.success === "true" || json.success === true)) {
        setStatus("success");
        f.reset();
      } else {
        setStatus("error");
        setError("The form server didn't accept the message. Try the email-client fallback below.");
      }
    } catch (err) {
      setStatus("error");
      setError("The form couldn't be reached from this browser (network or extension may be blocking it). Use the email-client fallback below.");
    }
  };

  const openInEmailClient = (e) => {
    e.preventDefault();
    const to = decodeTo();
    const body = "Name: \n\nAffiliation: \n\nType (Feedback / Question / Collaboration / Bug): \n\nMessage:\n";
    window.location.href = `mailto:${to}?subject=${encodeURIComponent(EMAIL_SUBJECT)}&body=${encodeURIComponent(body)}`;
  };

  const tryFallback = (e) => {
    e.preventDefault();
    fallbackSubmit();
  };

  return (
    <section className="feedback">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 07 / Contact · Feedback · Collaboration</span>
        </div>
        <h2 className="module-title">
          Found a bug, idea, or want to collaborate? Drop a note.
        </h2>
        <p className="module-sub">
          Messages are sent directly to the creator of the atlas with the subject
          &ldquo;CropGBWater – Feedback or Collaboration Request&rdquo;.
          Plain English is welcome. Your email is required; institutional affiliation is optional.
        </p>
      </div>

      <div className="feedback-grid">
        <form ref={formRef} className="feedback-form" onSubmit={handleSubmit} noValidate>
          <input type="hidden" name="_subject" value={EMAIL_SUBJECT} />
          <input type="hidden" name="_template" value="table" />
          <input type="hidden" name="_captcha" value="false" />
          <input type="hidden" name="_next"
                 value={(typeof window !== "undefined" ? window.location.origin + window.location.pathname : "") + "#contact"} />
          {/* Honeypot anti-spam */}
          <input type="text" name="_honey" style={{ display: "none" }} tabIndex={-1} autoComplete="off" />

          <label className="field">
            <span className="field-label">Your name</span>
            <input type="text" name="name" required maxLength={120} placeholder="Dr. Jane Doe" />
          </label>

          <label className="field">
            <span className="field-label">Your email <span style={{ color: "var(--accent)" }}>*</span></span>
            <input type="email" name="email" required maxLength={200} placeholder="jane@uni.edu" />
          </label>

          <label className="field">
            <span className="field-label">Affiliation (optional)</span>
            <input type="text" name="affiliation" maxLength={200} placeholder="University, lab, or organisation" />
          </label>

          <label className="field">
            <span className="field-label">Type of message</span>
            <select name="type" defaultValue="Feedback on the dashboard">
              <option>Feedback on the dashboard</option>
              <option>Methodological question</option>
              <option>Collaboration proposal</option>
              <option>Bug or data issue</option>
              <option>Press / media enquiry</option>
              <option>Other</option>
            </select>
          </label>

          <label className="field field-message">
            <span className="field-label">Message</span>
            <textarea name="message" required rows={7} maxLength={4000}
                      placeholder="What's on your mind? Be as specific as you'd like — paste links, mention dataset versions, or describe what you'd want the dashboard to do next." />
          </label>

          <div className="feedback-actions">
            <button type="submit" className="feedback-submit" disabled={status === "sending"}>
              {status === "sending" ? "Sending…" : "Send message"}
            </button>
            <a href="#" className="feedback-mailto" onClick={openInEmailClient}
               title="Opens your default email app with the subject pre-filled">
              Open in email client instead
            </a>
          </div>

          {status === "success" && (
            <div className="feedback-status feedback-success">
              <strong>Thanks — message sent.</strong> A reply will come back to you within a few days.
            </div>
          )}
          {status === "error" && (
            <div className="feedback-status feedback-error">
              <strong>Couldn&rsquo;t send via the form.</strong> {error}
              <div style={{ marginTop: 8, display: "flex", gap: 10, flexWrap: "wrap" }}>
                <a href="#" onClick={openInEmailClient}>Open in email client</a>
                <span style={{ color: "var(--ink40)" }}>·</span>
                <a href="#" onClick={tryFallback}>Retry via full page submit</a>
              </div>
            </div>
          )}
        </form>

        <aside className="feedback-aside">
          <div className="feedback-card">
            <h4>What we&rsquo;d love to hear</h4>
            <ul>
              <li>What works, what doesn&rsquo;t, what&rsquo;s missing</li>
              <li>Specific country / crop views you&rsquo;d want added</li>
              <li>Discrepancies between the dashboard and your own analysis</li>
              <li>Citation or attribution requests</li>
              <li>Collaboration on extending the model to new years or crops</li>
            </ul>
          </div>

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

          <div className="feedback-note">
            Privacy: messages are forwarded by FormSubmit.co to the project author and not stored on this site.
            The email address is not displayed on this page — only the subject line is visible to senders.
          </div>
        </aside>
      </div>
    </section>
  );
}

window.FeedbackPage = FeedbackPage;
