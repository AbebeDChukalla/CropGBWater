// feedback.jsx — Contact page wired to FormSubmit (https://formsubmit.co).
//
// FormSubmit posts the form fields directly to abedeme@gmail.com with the
// subject "CropGBWater – Feedback or Collaboration Request".
//
// First-time confirmation: when a form is first submitted with a new
// destination address, FormSubmit sends ONE activation email to that
// address. Click the link in that email; afterwards every submission is
// forwarded automatically — no signup, no API key.

const FORMSUBMIT_TARGET = "abedeme@gmail.com";
const FORM_ACTION = `https://formsubmit.co/${FORMSUBMIT_TARGET}`;
const EMAIL_SUBJECT = "CropGBWater – Feedback or Collaboration Request";

function FeedbackPage({ data }) {
  const [status, setStatus] = React.useState("idle"); // idle | sending | success | error
  const [error,  setError]  = React.useState(null);

  const handleSubmit = async (e) => {
    e.preventDefault();
    const f = e.currentTarget;
    setStatus("sending");
    setError(null);
    try {
      const fd = new FormData(f);
      // Ask FormSubmit to reply with JSON (so we can render a custom success state).
      const res = await fetch(`https://formsubmit.co/ajax/${FORMSUBMIT_TARGET}`, {
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
        setError(json.message || "Submission was rejected. Please try the mailto fallback below.");
      }
    } catch (err) {
      setStatus("error");
      setError(err.message || "Network error. Please try the mailto fallback below.");
    }
  };

  const mailtoFallback = () => {
    const body = encodeURIComponent(
      "Name: \n\nEmail: \n\nType (Feedback / Question / Collaboration / Bug): \n\nMessage:\n"
    );
    return `mailto:${FORMSUBMIT_TARGET}?subject=${encodeURIComponent(EMAIL_SUBJECT)}&body=${body}`;
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
          Messages go straight to <strong>{FORMSUBMIT_TARGET}</strong> with the subject
          &ldquo;CropGBWater – Feedback or Collaboration Request&rdquo;.
          Plain English is fine; institutional affiliations welcome but not required.
        </p>
      </div>

      <div className="feedback-grid">
        <form className="feedback-form" action={FORM_ACTION} method="POST" onSubmit={handleSubmit}>
          {/* FormSubmit config fields — all begin with underscore */}
          <input type="hidden" name="_subject" value={EMAIL_SUBJECT} />
          <input type="hidden" name="_template" value="table" />
          <input type="hidden" name="_captcha" value="false" />
          {/* Honeypot anti-spam — bots fill this, humans don't see it */}
          <input type="text" name="_honey" style={{ display: "none" }} tabIndex={-1} autoComplete="off" />

          <label className="field">
            <span className="field-label">Your name</span>
            <input type="text" name="name" required maxLength={120} placeholder="Dr. Jane Doe" />
          </label>

          <label className="field">
            <span className="field-label">Your email</span>
            <input type="email" name="email" required maxLength={200} placeholder="jane@uni.edu" />
          </label>

          <label className="field">
            <span className="field-label">Affiliation (optional)</span>
            <input type="text" name="affiliation" maxLength={200} placeholder="University, lab, or organisation" />
          </label>

          <label className="field">
            <span className="field-label">Type of message</span>
            <select name="type" defaultValue="Feedback">
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
            <a className="feedback-mailto" href={mailtoFallback()} title="Falls back to your email client if the form fails">
              Or open in email client
            </a>
          </div>

          {status === "success" && (
            <div className="feedback-status feedback-success">
              <strong>Thanks — message sent.</strong> A reply will come from {FORMSUBMIT_TARGET} (usually within a few days).
            </div>
          )}
          {status === "error" && (
            <div className="feedback-status feedback-error">
              <strong>Couldn&rsquo;t send via the form.</strong> {error}
              {" "}<a href={mailtoFallback()}>Use the mailto link instead</a>.
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
            <h4>Where else to reach us</h4>
            <ul className="feedback-links">
              <li><span className="action-key">DOI</span>
                  <a href={data?.meta?.paper_url} target="_blank" rel="noopener">{data?.meta?.paper_doi}</a></li>
              <li><span className="action-key">Data</span>
                  <a href={data?.meta?.data_url} target="_blank" rel="noopener">Zenodo · {data?.meta?.data_doi}</a></li>
              <li><span className="action-key">Code</span>
                  <a href="https://github.com/AbebeDChukalla/CropGBWater" target="_blank" rel="noopener">github.com/AbebeDChukalla/CropGBWater</a></li>
              <li><span className="action-key">Email</span>
                  <a href={mailtoFallback()}>{FORMSUBMIT_TARGET}</a></li>
            </ul>
          </div>

          <div className="feedback-note">
            Privacy: messages are forwarded by FormSubmit.co to a single mailbox and not stored on this site.
            If you&rsquo;d rather not use the web form, the mailto link opens your default email client with the
            subject pre-filled.
          </div>
        </aside>
      </div>
    </section>
  );
}

window.FeedbackPage = FeedbackPage;
