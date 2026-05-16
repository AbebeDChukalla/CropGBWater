// feedback.jsx — Contact page. Two destination addresses are displayed (and
// encoded), the message is delivered via the user's local email client (no
// third-party form host that can break). Includes simple inline validation
// and a loading state.

const TO_LIST_B64 = [
  "YWJlZGVtZUBnbWFpbC5jb20=",        // abedeme@gmail.com
  "YS5jaHVrYWxsYUB1bi1paGUub3Jn",    // a.chukalla@un-ihe.org
];
const decodeAll = () => TO_LIST_B64
  .map((s) => { try { return atob(s); } catch (_) { return ""; } })
  .filter(Boolean);
const EMAIL_SUBJECT = "CropGBWater – Feedback or Collaboration Request";

function FeedbackPage({ data }) {
  const [status, setStatus] = React.useState("idle"); // idle | sending | success
  const [errs, setErrs] = React.useState({});

  const validate = (fd) => {
    const e = {};
    const name    = (fd.get("name")    || "").toString().trim();
    const email   = (fd.get("email")   || "").toString().trim();
    const message = (fd.get("message") || "").toString().trim();
    if (!name) e.name = "Name is required.";
    if (!email) e.email = "Email is required.";
    else if (!/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(email)) e.email = "Please enter a valid email address.";
    if (!message || message.length < 10) e.message = "Message must be at least 10 characters.";
    return e;
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    const f = e.currentTarget;
    const fd = new FormData(f);
    const e2 = validate(fd);
    setErrs(e2);
    if (Object.keys(e2).length) return;

    setStatus("sending");

    const name        = fd.get("name").toString().trim();
    const email       = fd.get("email").toString().trim();
    const affiliation = (fd.get("affiliation") || "").toString().trim();
    const type        = (fd.get("type") || "Feedback").toString().trim();
    const message     = fd.get("message").toString().trim();

    const body = [
      `Name:        ${name}`,
      `Email:       ${email}`,
      affiliation ? `Affiliation: ${affiliation}` : null,
      `Type:        ${type}`,
      ``,
      `Message:`,
      message,
      ``,
      `— sent via the CropGBWater dashboard`,
    ].filter((l) => l !== null).join("\n");

    const recipients = decodeAll().join(",");
    const url = `mailto:${recipients}`
              + `?subject=${encodeURIComponent(EMAIL_SUBJECT)}`
              + `&body=${encodeURIComponent(body)}`;

    // Tiny delay so the loading state is visible (perception polish)
    setTimeout(() => {
      window.location.href = url;
      setStatus("success");
    }, 350);
  };

  return (
    <section className="feedback">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 07 / Contact · Feedback · Collaboration</span>
        </div>
        <h2 className="module-title">
          Reach the creator of the CropGBWater Atlas.
        </h2>
        <p className="module-sub">
          Messages go directly to the creator/author of the CropGBWater Atlas for
          feedback, scientific collaboration, data inquiries, educational use, or
          partnership opportunities. The subject line is set to
          &ldquo;{EMAIL_SUBJECT}&rdquo;. Your email is required; institutional
          affiliation is optional.
        </p>
      </div>

      <div className="feedback-grid">
        <form className="feedback-form" onSubmit={handleSubmit} noValidate>
          <label className="field">
            <span className="field-label">Your name <span className="req">*</span></span>
            <input type="text" name="name" required maxLength={120} placeholder="Dr. Jane Doe"
                   aria-invalid={!!errs.name}
                   className={errs.name ? "field-err" : ""} />
            {errs.name && <span className="field-msg">{errs.name}</span>}
          </label>

          <label className="field">
            <span className="field-label">Your email <span className="req">*</span></span>
            <input type="email" name="email" required maxLength={200} placeholder="jane@uni.edu"
                   aria-invalid={!!errs.email}
                   className={errs.email ? "field-err" : ""} />
            {errs.email && <span className="field-msg">{errs.email}</span>}
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
              <option>Data inquiry</option>
              <option>Educational / teaching use</option>
              <option>Partnership opportunity</option>
              <option>Bug or data issue</option>
              <option>Press / media enquiry</option>
              <option>Other</option>
            </select>
          </label>

          <label className="field field-message">
            <span className="field-label">Message <span className="req">*</span></span>
            <textarea name="message" required rows={7} maxLength={4000}
                      aria-invalid={!!errs.message}
                      className={errs.message ? "field-err" : ""}
                      placeholder="What's on your mind? Be as specific as you'd like — paste links, mention dataset versions, or describe what you'd want the dashboard to do next." />
            {errs.message && <span className="field-msg">{errs.message}</span>}
            {/* Honeypot anti-spam */}
            <input type="text" name="_honey" style={{ display: "none" }} tabIndex={-1} autoComplete="off" />
          </label>

          <div className="feedback-actions">
            <button type="submit" className="feedback-submit" disabled={status === "sending"}>
              {status === "sending" ? (
                <span className="spinner" aria-hidden="true" />
              ) : null}
              {status === "sending" ? "Opening email…" : "Send message"}
            </button>
            <span className="feedback-mailto-hint">Opens in your default email app</span>
          </div>

          {status === "success" && (
            <div className="feedback-status feedback-success">
              <strong>Email client opened.</strong> Click &ldquo;Send&rdquo; in your email
              app to deliver the message. If nothing opened, check your browser&rsquo;s pop-up
              settings or write directly to one of the addresses on the right.
            </div>
          )}
        </form>

        <aside className="feedback-aside">
          <div className="feedback-card">
            <h4>Contact addresses</h4>
            <ul className="feedback-links">
              <li><span className="action-key">Author</span>
                  <button type="button" className="reveal-email" data-email-b64={TO_LIST_B64[0]}
                          onClick={revealAndCopy}>abedeme@gmail.com</button></li>
              <li><span className="action-key">Inst.</span>
                  <button type="button" className="reveal-email" data-email-b64={TO_LIST_B64[1]}
                          onClick={revealAndCopy}>a.chukalla@un-ihe.org</button></li>
            </ul>
            <span className="feedback-note-inline">Click to copy to clipboard.</span>
          </div>

          <div className="feedback-card">
            <h4>What we&rsquo;d love to hear</h4>
            <ul>
              <li>What works, what doesn&rsquo;t, what&rsquo;s missing</li>
              <li>Specific country / crop views you&rsquo;d want added</li>
              <li>Discrepancies between the dashboard and your own analysis</li>
              <li>Citation or attribution requests</li>
              <li>Collaboration on extending the model to new years or crops</li>
              <li>Educational use in courses or workshops</li>
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
            Privacy: submissions are handed to your local email app and not stored on this
            site or sent through any third-party form service. Anti-spam: a honeypot field
            blocks automated submissions.
          </div>
        </aside>
      </div>
    </section>
  );
}

// Reveal + copy-to-clipboard for the contact email buttons
function revealAndCopy(e) {
  const b64 = e.currentTarget.getAttribute("data-email-b64") || "";
  let addr = "";
  try { addr = atob(b64); } catch (_) {}
  if (!addr) return;
  try {
    navigator.clipboard.writeText(addr);
    const orig = e.currentTarget.textContent;
    e.currentTarget.textContent = "Copied ✓";
    setTimeout(() => { e.currentTarget.textContent = orig; }, 1500);
  } catch (_) {
    // Fallback: open mailto if clipboard not available
    window.location.href = `mailto:${addr}`;
  }
}

window.FeedbackPage = FeedbackPage;
