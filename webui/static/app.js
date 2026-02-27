const API_BASE = "/api";
const POLL_INTERVAL = 2000;
const MAX_POLL_TIME = 600000;

const form = document.getElementById("analysisForm");
const submitBtn = document.getElementById("submitBtn");
const quickTestBtn = document.getElementById("quickTestBtn");
const loadJobBtn = document.getElementById("loadJobBtn");
const existingJobIdInput = document.getElementById("existingJobId");
const recentJobsSelect = document.getElementById("recentJobsSelect");
const backBtn = document.getElementById("backBtn");
const eventTypeSelect = document.getElementById("eventType");
const genomeSelect = document.getElementById("genome");
const inputTypeRadios = document.querySelectorAll('input[name="input_type"]');
const rmatsInput = document.getElementById("rmatsInput");
const misoInput = document.getElementById("misoInput");
const coordinatesInput = document.getElementById("coordinatesInput");
const resultsDiv = document.getElementById("results");
const liveLogsEl = document.getElementById("liveLogs");
const liveLogsWrapEl = document.getElementById("liveLogsWrap");

let currentJobId = null;
let pollStartTime = null;
let currentOutputDir = null;
let suppressHistoryPush = false;
let loadJobInFlight = false;

document.addEventListener("DOMContentLoaded", () => {
  loadGenomes();
  loadRecentJobs();
  setupEventListeners();
  loadFromStorage();
  applyLocationState();
  window.addEventListener("popstate", () => {
    applyLocationState();
  });
});

async function loadRecentJobs() {
  if (!recentJobsSelect) return;
  try {
    const response = await fetch(`${API_BASE}/jobs?limit=30`);
    const data = await response.json();
    if (!response.ok || !data.success) {
      recentJobsSelect.innerHTML = '<option value="">-- Recent Jobs (unavailable) --</option>';
      return;
    }
    recentJobsSelect.innerHTML = '<option value="">-- Recent Jobs --</option>';
    const jobs = data.jobs || [];
    for (const job of jobs) {
      const option = document.createElement("option");
      option.value = job.job_id;
      option.textContent = `${job.job_id} (${job.name})`;
      recentJobsSelect.appendChild(option);
    }
    if (jobs.length === 0) {
      recentJobsSelect.innerHTML = '<option value="">-- No Jobs Found --</option>';
    }
  } catch (error) {
    console.error("Error loading recent jobs:", error);
    recentJobsSelect.innerHTML = '<option value="">-- Recent Jobs (error) --</option>';
  }
}

async function loadGenomes() {
  try {
    const response = await fetch(`${API_BASE}/genomes`);
    const data = await response.json();
    const genomes = data.genomes || [];
    if (!Array.isArray(genomes) || genomes.length === 0) {
      return;
    }

    genomeSelect.innerHTML = "";
    let selected = false;
    for (const genome of genomes) {
      const option = document.createElement("option");
      option.value = genome.code;
      const availability = genome.available ? "" : " (FASTA missing)";
      option.textContent = `${genome.code} - ${genome.organism}${availability}`;
      if (!selected && (genome.default || genome.code === "hg19")) {
        option.selected = true;
        selected = true;
      }
      genomeSelect.appendChild(option);
    }
    if (!selected && genomeSelect.options.length > 0) {
      genomeSelect.options[0].selected = true;
    }
  } catch (error) {
    console.error("Error loading genomes:", error);
    showError("Failed to load genomes. Please refresh the page.");
  }
}

function setupEventListeners() {
  inputTypeRadios.forEach((radio) => {
    radio.addEventListener("change", (e) => {
      if (e.target.value === "rmats") {
        rmatsInput.style.display = "block";
        misoInput.style.display = "none";
        coordinatesInput.style.display = "none";
        document.getElementById("rmatsFile").required = true;
        document.getElementById("misoFile").required = false;
        document.getElementById("upFile").required = false;
        document.getElementById("dnFile").required = false;
        document.getElementById("bgFile").required = false;
      } else if (e.target.value === "miso") {
        rmatsInput.style.display = "none";
        misoInput.style.display = "block";
        coordinatesInput.style.display = "none";
        document.getElementById("rmatsFile").required = false;
        document.getElementById("misoFile").required = true;
        document.getElementById("upFile").required = false;
        document.getElementById("dnFile").required = false;
        document.getElementById("bgFile").required = false;
      } else {
        rmatsInput.style.display = "none";
        misoInput.style.display = "none";
        coordinatesInput.style.display = "block";
        document.getElementById("rmatsFile").required = false;
        document.getElementById("misoFile").required = false;
        document.getElementById("upFile").required = true;
        document.getElementById("dnFile").required = true;
        document.getElementById("bgFile").required = true;
      }
    });
  });

  form.addEventListener("submit", handleFormSubmit);
  if (quickTestBtn) {
    quickTestBtn.addEventListener("click", runQuickTest);
  }
  if (loadJobBtn) {
    loadJobBtn.addEventListener("click", loadExistingJob);
  }
  if (recentJobsSelect && existingJobIdInput) {
    recentJobsSelect.addEventListener("change", () => {
      const value = recentJobsSelect.value || "";
      if (value) {
        existingJobIdInput.value = value;
      }
    });
  }
  if (backBtn) {
    backBtn.addEventListener("click", resetUI);
  }
  form.addEventListener("change", saveToStorage);
}

async function handleFormSubmit(e) {
  e.preventDefault();
  if (!form.checkValidity()) {
    form.reportValidity();
    return;
  }

  const formData = new FormData(form);
  const inputType = document.querySelector('input[name="input_type"]:checked').value;

  if (inputType === "rmats" && (!formData.get("rmats_file") || formData.get("rmats_file").size === 0)) {
    showError("Please upload an rMATS file");
    return;
  }
  if (inputType === "miso" && (!formData.get("miso_file") || formData.get("miso_file").size === 0)) {
    showError("Please upload a MISO file");
    return;
  }
  if (inputType === "coordinates") {
    if (!formData.get("up_file") || formData.get("up_file").size === 0) {
      showError("Please upload upregulated exons file");
      return;
    }
    if (!formData.get("down_file") || formData.get("down_file").size === 0) {
      showError("Please upload downregulated exons file");
      return;
    }
    if (!formData.get("bg_file") || formData.get("bg_file").size === 0) {
      showError("Please upload background exons file");
      return;
    }
  }
  if (!formData.get("known_motifs_file") || formData.get("known_motifs_file").size === 0) {
    showError("Please upload a known motifs file");
    return;
  }

  setRunButtonsDisabled(true, "Submitting...");
  try {
    const response = await fetch(`${API_BASE}/submit`, { method: "POST", body: formData });
    const data = await response.json();
    if (!response.ok || !data.success) {
      showError(data.error || "Failed to submit analysis");
      setRunButtonsDisabled(false);
      return;
    }

    currentJobId = data.job_id;
    currentOutputDir = data.output_dir || null;
    pollStartTime = Date.now();
    updateJobId(currentJobId);
    setUrlForJob(currentJobId);
    showMessage(`Analysis submitted. Job ID: ${currentJobId}`, "info");
    pollJobStatus();
    form.style.display = "none";
  } catch (error) {
    console.error("Error submitting form:", error);
    showError(`Error submitting analysis: ${error.message}`);
    setRunButtonsDisabled(false);
  }
}

async function runQuickTest() {
  setRunButtonsDisabled(true, "Starting Test...");
  try {
    const eventType = eventTypeSelect.value || "se";
    const genome = genomeSelect.value || "hg19";
    const response = await fetch(`${API_BASE}/quick-test/run`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ event_type: eventType, genome }),
    });
    const data = await response.json();
    if (!response.ok || !data.success) {
      if (data.missing_files && data.quick_test_dir) {
        const missingText = data.missing_files.map((p) => `- ${p}`).join("\n");
        showError(`Quick test files missing.\nExpected directory: ${data.quick_test_dir}\n${missingText}`);
      } else {
        showError(data.error || "Failed to start one-click test");
      }
      setRunButtonsDisabled(false);
      return;
    }

    currentJobId = data.job_id;
    currentOutputDir = data.output_dir || null;
    pollStartTime = Date.now();
    updateJobId(currentJobId);
    setUrlForJob(currentJobId);
    showMessage(`One-click test submitted. Job ID: ${currentJobId}`, "info");
    pollJobStatus();
    form.style.display = "none";
  } catch (error) {
    console.error("Error starting one-click test:", error);
    showError(`Error starting one-click test: ${error.message}`);
    setRunButtonsDisabled(false);
  }
}

async function loadExistingJob(forcedJobId = null, skipHistory = false) {
  if (loadJobInFlight) {
    return;
  }
  const picked = recentJobsSelect && recentJobsSelect.value ? recentJobsSelect.value : "";
  const raw = forcedJobId || picked || ((existingJobIdInput && existingJobIdInput.value) ? existingJobIdInput.value : "");
  const jobId = String(raw).trim();
  if (!jobId) {
    showError("Enter a job ID first.");
    return;
  }

  loadJobInFlight = true;
  if (loadJobBtn) {
    loadJobBtn.disabled = true;
    loadJobBtn.textContent = "Loading...";
  }
  setRunButtonsDisabled(true, "Loading...");
  try {
    const selectedText = recentJobsSelect && recentJobsSelect.selectedOptions && recentJobsSelect.selectedOptions[0]
      ? (recentJobsSelect.selectedOptions[0].textContent || "")
      : "";
    const match = selectedText.match(/\(([^)]+)\)/);
    const selectedName = match ? match[1] : "";
    const alt = jobId.startsWith("quicktest_") ? jobId.replace(/^quicktest_/, "") : `quicktest_${jobId}`;
    const idCandidates = [];
    [jobId, alt, selectedName].forEach((x) => {
      if (x && !idCandidates.includes(x)) {
        idCandidates.push(x);
      }
    });

    let statusResp = null;
    let statusData = null;
    let resolvedJobId = jobId;
    for (const candidate of idCandidates) {
      const resp = await fetch(`${API_BASE}/status/${candidate}`);
      const data = await resp.json();
      if (resp.ok && data.success) {
        statusResp = resp;
        statusData = data;
        resolvedJobId = candidate;
        break;
      }
      if (!statusResp) {
        statusResp = resp;
        statusData = data;
      }
    }

    if (!statusResp.ok || !statusData.success) {
      showError(statusData.error || `Job ${jobId} not found`);
      setRunButtonsDisabled(false);
      return;
    }

    currentJobId = resolvedJobId;
    currentOutputDir = statusData.output_dir || null;
    pollStartTime = Date.now();
    updateJobId(currentJobId);
    if (!skipHistory) {
      setUrlForJob(currentJobId);
    }
    resultsDiv.style.display = "block";
    form.style.display = "none";
    await refreshJobLogs();
    updateStatusDisplay(statusData.status, statusData.message || "");

    if (statusData.status === "completed") {
      const resultsResp = await fetch(`${API_BASE}/results/${currentJobId}`);
      const resultsData = await resultsResp.json();
      if (!resultsResp.ok || !resultsData.success) {
        await showResults([], currentOutputDir, resultsData.error || "Completed. Summary table unavailable.");
      } else {
        await showResults(resultsData.results || [], currentOutputDir);
      }
      setRunButtonsDisabled(false);
      return;
    }

    if (statusData.status === "failed") {
      showError(statusData.message || "This job failed.");
      setRunButtonsDisabled(false);
      return;
    }

    pollJobStatus();
  } catch (error) {
    console.error("Error loading existing job:", error);
    showError(`Error loading job: ${error.message}`);
    setRunButtonsDisabled(false);
  } finally {
    loadJobInFlight = false;
    if (loadJobBtn) {
      loadJobBtn.disabled = false;
      loadJobBtn.textContent = "Load Existing Job";
    }
  }
}

async function pollJobStatus() {
  if (!currentJobId) return;
  if (Date.now() - pollStartTime > MAX_POLL_TIME) {
    showError("Analysis timeout (10 minutes).");
    resetUI();
    return;
  }

  try {
    const response = await fetch(`${API_BASE}/status/${currentJobId}`);
    const data = await response.json();
    await refreshJobLogs();

    if (!data.success) {
      showError("Error checking job status");
      resetUI();
      return;
    }

    currentOutputDir = data.output_dir || currentOutputDir;
    updateJobId(currentJobId);
    updateStatusDisplay(data.status, data.message);
    if (data.status === "completed") {
      const resultsResponse = await fetch(`${API_BASE}/results/${currentJobId}`);
      const resultsData = await resultsResponse.json();
      if (!resultsResponse.ok || !resultsData.success) {
        await showResults([], currentOutputDir, resultsData.error || "Completed. Summary table unavailable.");
      } else {
        await showResults(resultsData.results, currentOutputDir);
      }
      setRunButtonsDisabled(false);
    } else if (data.status === "failed") {
      showError(`Analysis failed: ${data.message}`);
      updateStatusDisplay("failed", data.message || "Pipeline failed");
      setRunButtonsDisabled(false);
    } else {
      setTimeout(pollJobStatus, POLL_INTERVAL);
    }
  } catch (error) {
    console.error("Error polling job status:", error);
    showError(`Error checking job status: ${error.message}`);
    setRunButtonsDisabled(false);
  }
}

async function refreshJobLogs() {
  if (!currentJobId || !liveLogsEl) return;
  try {
    const response = await fetch(`${API_BASE}/logs/${currentJobId}?tail=140`);
    const data = await response.json();
    if (!data.success) return;

    if (liveLogsWrapEl) {
      liveLogsWrapEl.style.display = "block";
    }
    const text = (data.lines || []).join("\n");
    const prior = liveLogsEl.textContent || "";
    liveLogsEl.textContent = text || "Waiting for logs...";
    if (text && text !== prior) {
      liveLogsEl.scrollTop = liveLogsEl.scrollHeight;
    }
  } catch (error) {
    console.error("Error fetching logs:", error);
  }
}

function updateStatusDisplay(status, message) {
  const statusSpan = document.querySelector("#results #status");
  if (!statusSpan) return;

  const badge = document.createElement("span");
  badge.className = `status-badge ${status}`;
  badge.textContent = status.charAt(0).toUpperCase() + status.slice(1);
  statusSpan.innerHTML = "";
  statusSpan.appendChild(badge);
  if (message) {
    statusSpan.appendChild(document.createTextNode(` - ${message}`));
  }
}

function updateJobId(jobId) {
  const el = document.getElementById("jobId");
  if (el) {
    el.textContent = jobId || "";
  }
}

async function showResults(results, outputDir, note = "") {
  resultsDiv.style.display = "block";
  if (liveLogsWrapEl) {
    liveLogsWrapEl.style.display = "block";
  }
  const content = document.getElementById("resultsContent");
  let html = "";
  if (note) {
    html += `<div class="info">${note}</div>`;
  }

  if (!results || results.length === 0) {
    html += "<p>No significant motifs found or summary table unavailable.</p>";
  } else {
    const topResults = results.slice(0, 10);
    html += "<h3>Top Motifs</h3>";
    html += "<table><thead><tr><th>Motif</th><th>p(up vs bg)</th><th>-log10 p(up)</th><th>p(dn vs bg)</th><th>-log10 p(dn)</th></tr></thead><tbody>";
    for (const result of topResults) {
      html += `<tr><td><strong>${result.rbp}</strong></td><td>${parseFloat(result.pval_up_vs_bg).toExponential(2)}</td><td>${parseFloat(result.log_pval_up).toFixed(2)}</td><td>${parseFloat(result.pval_dn_vs_bg).toExponential(2)}</td><td>${parseFloat(result.log_pval_dn).toFixed(2)}</td></tr>`;
    }
    html += "</tbody></table>";
    if (results.length > topResults.length) {
      html += `<p>Showing top ${topResults.length} of ${results.length} motifs.</p>`;
    }
  }
  html += `<div style="margin-top: 18px;"><h3>Result Folder</h3><pre>${outputDir || "(unknown path)"}</pre><p>Logs: ${(outputDir || "(unknown path)") + "/log.motifMap.txt"}</p></div>`;

  content.innerHTML = html;
  updateStatusDisplay("completed", "Analysis complete");
}

function showError(message) {
  const existing = resultsDiv.querySelectorAll(".error");
  existing.forEach((el) => el.remove());
  const div = document.createElement("div");
  div.className = "error";
  div.textContent = message;
  resultsDiv.style.display = "block";
  resultsDiv.insertBefore(div, resultsDiv.firstChild);
  setTimeout(() => div.remove(), 7000);
}

function showMessage(message, type = "info") {
  const div = document.createElement("div");
  div.className = type;
  div.textContent = message;
  resultsDiv.style.display = "block";
  resultsDiv.insertBefore(div, resultsDiv.firstChild);
}

function resetUI() {
  currentJobId = null;
  currentOutputDir = null;
  setRunButtonsDisabled(false);
  updateJobId("");
  setUrlForHome();
  form.style.display = "block";
  resultsDiv.style.display = "none";
  document.getElementById("resultsContent").innerHTML = "";
  if (liveLogsEl) {
    liveLogsEl.textContent = "";
  }
  loadRecentJobs();
}

function setRunButtonsDisabled(disabled, submitLabel = "Run Analysis") {
  submitBtn.disabled = disabled;
  submitBtn.textContent = disabled ? submitLabel : "Run Analysis";
  if (quickTestBtn) {
    quickTestBtn.disabled = disabled;
    quickTestBtn.textContent = disabled ? "Please wait..." : "Run One-Click Test";
  }
}

function saveToStorage() {
  const selectedInput = document.querySelector('input[name="input_type"]:checked');
  const formState = {
    eventType: eventTypeSelect.value,
    genome: genomeSelect.value,
    inputType: selectedInput ? selectedInput.value : "rmats",
    label: document.getElementById("label").value,
  };
  localStorage.setItem("rmapsFormState", JSON.stringify(formState));
}

function loadFromStorage() {
  const saved = localStorage.getItem("rmapsFormState");
  if (!saved) return;
  try {
    const state = JSON.parse(saved);
    if (state.eventType) eventTypeSelect.value = state.eventType;
    if (state.inputType) {
      const radio = document.querySelector(`input[name="input_type"][value="${state.inputType}"]`);
      if (radio) {
        radio.checked = true;
        radio.dispatchEvent(new Event("change"));
      }
    }
    if (state.label) document.getElementById("label").value = state.label;
  } catch (error) {
    console.error("Error loading saved form state:", error);
  }
}

function setUrlForHome() {
  if (suppressHistoryPush) return;
  const url = new URL(window.location.href);
  url.searchParams.delete("job");
  window.history.pushState({ view: "home" }, "", url.toString());
}

function setUrlForJob(jobId) {
  if (!jobId || suppressHistoryPush) return;
  const url = new URL(window.location.href);
  url.searchParams.set("job", jobId);
  window.history.pushState({ view: "job", jobId }, "", url.toString());
}

function applyLocationState() {
  const url = new URL(window.location.href);
  const jobId = (url.searchParams.get("job") || "").trim();
  suppressHistoryPush = true;
  try {
    if (jobId) {
      if (existingJobIdInput) {
        existingJobIdInput.value = jobId;
      }
      loadExistingJob(jobId, true);
    } else {
      currentJobId = null;
      currentOutputDir = null;
      setRunButtonsDisabled(false);
      updateJobId("");
      form.style.display = "block";
      resultsDiv.style.display = "none";
      document.getElementById("resultsContent").innerHTML = "";
      if (liveLogsEl) {
        liveLogsEl.textContent = "";
      }
    }
  } finally {
    suppressHistoryPush = false;
  }
}
