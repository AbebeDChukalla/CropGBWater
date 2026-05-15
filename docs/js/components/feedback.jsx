// feedback.jsx — Contact page. The destination email is hidden from the page:
// base64-encoded in JS, decoded only when the user clicks Send.
//
// Delivery is via mailto: — when the user submits, we read all the form
// fields, build a structured plain-text body, and hand it off to the user's
// default email client. This works regardless of network state, CORS, or
// browser extensions — no third-party form service to depend on.

const TO_B64 = "YWJlZGVtZUBnbWFpbC5jb20=";       // base64 of the destination
const decodeTo = () => { try { return atob(TO_B64); } catch (_) { return ""; } };
const EMAIL_SUBJECT = "CropGBWater – Feedback or Collaboration Request";

function FeedbackPage({ data }) {
  const [sent, setSent] = React.useState(false);

  const handleSubmit = (e) => {
    e.preventDefault();
    const f = e.currentTarget;
    const fd = new FormData(f);
    const name        = (fd.get("name")        || "").toString().trim();
    const email       = (fd.get("email")       || "").toString().trim();
    const affiliation = (fd.get("affiliation") || "").toString().trim();
    const type        = (fd.get("type")        || "Feedback").toString().trim();
    const message     = (fd.get("message")     || "").toString().trim();

    if (!name || !email || !message) return;

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

    const url = `mailto:${decodeTo()}`
              + `?subject=${encodeURIComponent(EMAIL_SUBJECT)}`
              + `&body=${encodeURIComponent(body)}`;

    // Open the user's default email client with everything pre-filled
    window.location.href = url;
    // Mark sent UI on the page (user still needs to click Send in their mail app)
    setSent(true);
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
        <form className="feedback-form" onSubmit={handleSubmit} noValidate>
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
            <button type="submit" className="feedback-submit">Send message</button>
            <span className="feedback-mailto-hint">Opens in your default email app</span>
          </div>

          {sent && (
            <div className="feedback-status feedback-success">
              <strong>Email client opened.</strong> Click &ldquo;Send&rdquo; in your email app to deliver the message.
              {" "}If nothing opened, check your browser&rsquo;s pop-up settings or use your normal email program directly.
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
            Privacy: the destination email is not displayed on this page.
            Submissions are handed to your local email app and not stored on this site
            or sent through any third-party form service.
          </div>
        </aside>
      </div>
    </section>
  );
}

window.FeedbackPage = FeedbackPage;
