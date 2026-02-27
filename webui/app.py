from __future__ import annotations

import logging
import math
import os
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from threading import Lock
from typing import Any

from flask import Flask, jsonify, render_template, request

from motif_map_core import EVENT_SPECS, run_motif_map


LOGGER = logging.getLogger(__name__)


SUPPORTED_GENOMES = {
    "hg38": {"organism": "Homo sapiens", "default": True},
    "hg19": {"organism": "Homo sapiens", "default": False},
    "mm10": {"organism": "Mus musculus", "default": False},
    "bosTau9": {"organism": "Bos taurus", "default": False},
    "dm6": {"organism": "Drosophila melanogaster", "default": False},
    "dm3": {"organism": "Drosophila melanogaster", "default": False},
    "rn6": {"organism": "Rattus norvegicus", "default": False},
    "ce11": {"organism": "Caenorhabditis elegans", "default": False},
    "danRer10": {"organism": "Danio rerio", "default": False},
    "danRer11": {"organism": "Danio rerio", "default": False},
    "galGal5": {"organism": "Gallus gallus", "default": False},
    "galGal6": {"organism": "Gallus gallus", "default": False},
    "araTha1": {"organism": "Arabidopsis thaliana", "default": False},
    "oSa7": {"organism": "Oryza sativa", "default": False},
    "susScr11": {"organism": "Sus scrofa", "default": False},
    "xenLae2": {"organism": "Xenopus laevis", "default": False},
    "xenTro7": {"organism": "Xenopus tropicalis", "default": False},
    "xenTro9": {"organism": "Xenopus tropicalis", "default": False},
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _tail_lines(path: Path, max_lines: int) -> list[str]:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            return handle.readlines()[-max_lines:]
    except Exception:
        return []


def _important_log_lines(lines: list[str]) -> list[str]:
    levels = ("| INFO", "| WARNING", "| ERROR", "| CRITICAL")
    drop_markers = (
        '"GET /api/',
        '"POST /api/',
        '"GET /static/',
        'HTTP/1.1" 200',
        'HTTP/1.1" 304',
        'HTTP/1.1" 400',
        'HTTP/1.1" 404',
        'HTTP/1.1" 500',
    )
    keep: list[str] = []
    suppressed = {
        "plotregions": 0,
        "pyx_locator": 0,
        "native_png_failed": 0,
    }
    rendered_maps = 0
    last_sep = False
    for raw in lines:
        line = raw.rstrip()
        if not line:
            continue
        if any(marker in line for marker in drop_markers):
            continue
        if "Done plotRegions function" in line:
            suppressed["plotregions"] += 1
            continue
        if (
            "PyX executes kpsewhich" in line
            or "PyX filelocator" in line
            or "PyX executes tex" in line
            or "Evaluating lazy import" in line
        ):
            suppressed["pyx_locator"] += 1
            continue
        if "Native PNG export failed" in line:
            suppressed["native_png_failed"] += 1
            # Keep only first instance for context.
            if suppressed["native_png_failed"] > 1:
                continue
        if "Done drawing acutal plots" in line:
            rendered_maps += 1
            if rendered_maps % 10 == 0:
                keep.append(f"[progress] rendered motif maps: {rendered_maps}")
            continue
        if "================================" in line:
            if last_sep:
                continue
            last_sep = True
            keep.append(line)
            continue
        last_sep = False
        # Keep explicit leveled logs and classic legacy timestamp lines.
        if any(level in line for level in levels) or (len(line) > 4 and line[:4].isdigit()):
            keep.append(line)
            continue
        # Fallback: keep any non-empty non-noise line so live output does not look frozen.
        keep.append(line)
    if rendered_maps > 0 and rendered_maps % 10 != 0:
        keep.append(f"[progress] rendered motif maps: {rendered_maps}")
    if suppressed["plotregions"] > 0:
        keep.append(f"[summary] suppressed {suppressed['plotregions']} repeated plot-region lines")
    if suppressed["pyx_locator"] > 0:
        keep.append(f"[summary] suppressed {suppressed['pyx_locator']} repeated PyX locator lines")
    if suppressed["native_png_failed"] > 1:
        keep.append(f"[summary] suppressed {suppressed['native_png_failed'] - 1} repeated PNG-export warnings")
    return keep


def _normalize_event(value: str) -> str:
    return value.strip().lower()


def _create_job_dir(base: Path, prefix: str = "") -> Path:
    job_id = str(uuid.uuid4())[:8]
    name = f"{prefix}{job_id}" if prefix else job_id
    path = base / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def create_app() -> Flask:
    app = Flask(
        __name__,
        template_folder=str(Path(__file__).resolve().parent / "templates"),
        static_folder=str(Path(__file__).resolve().parent / "static"),
        static_url_path="/static",
    )
    app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024

    repo = _repo_root()
    app.config["REPO_ROOT"] = repo
    app.config["FASTA_ROOT"] = Path(os.environ.get("RMAPS_FASTA_ROOT", str(repo / "genomedata")))
    app.config["RESULTS_FOLDER"] = Path(os.environ.get("RMAPS_RESULTS_DIR", str(repo / "results")))
    app.config["QUICK_TEST_DIR"] = Path(os.environ.get("RMAPS_QUICKTEST_DIR", str(repo / "testData")))
    app.config["RESULTS_FOLDER"].mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
    app.jobs: dict[str, dict[str, Any]] = {}
    app.jobs_lock = Lock()
    app.executor = ThreadPoolExecutor(max_workers=2, thread_name_prefix="rmaps-web")

    register_routes(app)
    return app


def register_routes(app: Flask) -> None:
    def resolve_job(job_id: str) -> dict[str, Any] | None:
        with app.jobs_lock:
            in_memory = app.jobs.get(job_id)
        if in_memory:
            return dict(in_memory)

        results_root = Path(app.config["RESULTS_FOLDER"])
        candidates = [results_root / job_id, results_root / f"quicktest_{job_id}"]
        for path in candidates:
            if path.exists() and path.is_dir():
                log_exists = (path / "log.motifMap.txt").exists()
                status = "completed" if log_exists else "unknown"
                return {
                    "status": status,
                    "event_type": "unknown",
                    "genome": "unknown",
                    "created": datetime.fromtimestamp(path.stat().st_mtime),
                    "output_dir": str(path),
                    "result_dir": str(path),
                    "error": None,
                }
        return None

    def quick_test_files() -> dict[str, Path]:
        root = Path(app.config["QUICK_TEST_DIR"])
        rmats_candidates = [
            root / "motifMap.testEvents.rMATS.txt",
            root / "SE.MATS.ReadsOnTargetAndJunctionCounts.txt",
            root / "motifMap.testEvents.txt",
        ]
        known_candidates = [
            root / "knownMotifs.human.mouse.txt",
            Path(app.config["REPO_ROOT"]) / "data" / "knownMotifs.human.mouse.txt",
            Path(app.config["REPO_ROOT"]) / "data" / "short.knownMotifs.human.mouse.txt",
        ]
        motifs_candidates = [
            root / "ESRP.like.motif.txt",
            Path(app.config["REPO_ROOT"]) / "data" / "ESRP.like.motif.txt",
        ]

        def first_existing(paths: list[Path]) -> Path:
            for p in paths:
                if p.exists() and p.is_file():
                    return p
            return paths[0]

        return {
            "root": root,
            "rmats": first_existing(rmats_candidates),
            "known_motifs": first_existing(known_candidates),
            "motifs": first_existing(motifs_candidates),
        }

    def missing_quick_test_paths() -> list[str]:
        q = quick_test_files()
        out = []
        for key in ("rmats", "known_motifs", "motifs"):
            p = q[key]
            if not p.exists() or not p.is_file():
                out.append(str(p))
        return out

    def run_job(job_id: str) -> None:
        with app.jobs_lock:
            job = app.jobs.get(job_id)
            if not job:
                return
            job["status"] = "running"
            job = dict(job)

        params = job["params"]
        input_files = job["input_files"]
        event = _normalize_event(job["event_type"])

        code = run_motif_map(
            event=event,
            known_motifs=Path(input_files["known_motifs"]),
            motifs=input_files.get("motifs", "NA"),
            fasta_root=Path(app.config["FASTA_ROOT"]),
            genome=job["genome"],
            output=Path(job["output_dir"]),
            rmats=input_files.get("rmats", "NA"),
            miso=input_files.get("miso", "NA"),
            up=input_files.get("up", "NA"),
            down=input_files.get("down", "NA"),
            background=input_files.get("bg", "NA"),
            label=params["rbp_label"],
            intron=params["intron_len"],
            exon=params["exon_len"],
            window=params["window_size"],
            step=params["step_size"],
            sig_fdr=params["sig_fdr"],
            sig_delta_psi=params["sig_delta_psi"],
            separate=False,
        )

        with app.jobs_lock:
            if job_id not in app.jobs:
                return
            if code == 0:
                app.jobs[job_id]["status"] = "completed"
                app.jobs[job_id]["result_dir"] = job["output_dir"]
            else:
                app.jobs[job_id]["status"] = "failed"
                app.jobs[job_id]["error"] = f"Pipeline failed with exit code {code}"

    @app.get("/")
    def index() -> str:
        return render_template("index.html")

    @app.get("/api/events")
    def get_events():
        events = []
        for code in EVENT_SPECS:
            events.append(
                {
                    "code": code,
                    "name": EVENT_SPECS[code].name,
                    "description": f"{EVENT_SPECS[code].name} event",
                }
            )
        return jsonify({"success": True, "events": events}), 200

    @app.get("/api/genomes")
    def get_genomes():
        genomes = []
        fasta_root: Path = app.config["FASTA_ROOT"]
        available = set()
        if fasta_root.exists():
            for child in fasta_root.iterdir():
                if child.is_dir():
                    fa = child / f"{child.name}.fa"
                    if fa.exists():
                        available.add(child.name)

        keys = set(SUPPORTED_GENOMES.keys()) | available
        for code in sorted(keys):
            meta = SUPPORTED_GENOMES.get(code, {"organism": "Unknown", "default": False})
            genomes.append(
                {
                    "code": code,
                    "name": code,
                    "organism": meta["organism"],
                    "default": meta["default"],
                    "available": code in available,
                }
            )
        return jsonify({"success": True, "genomes": genomes}), 200

    @app.get("/api/quick-test/config")
    def quick_test_config():
        q = quick_test_files()
        missing = missing_quick_test_paths()
        return (
            jsonify(
                {
                    "success": True,
                    "quick_test_dir": str(q["root"]),
                    "required_files": {
                        "rmats": str(q["rmats"]),
                        "known_motifs": str(q["known_motifs"]),
                        "motifs": str(q["motifs"]),
                    },
                    "missing_files": missing,
                    "ready": len(missing) == 0,
                }
            ),
            200,
        )

    @app.post("/api/quick-test/run")
    def run_quick_test():
        missing = missing_quick_test_paths()
        q = quick_test_files()
        if missing:
            return (
                jsonify(
                    {
                        "success": False,
                        "error": "Quick test files are missing",
                        "quick_test_dir": str(q["root"]),
                        "missing_files": missing,
                    }
                ),
                400,
            )

        payload = request.get_json(silent=True) or {}
        event_type = _normalize_event(payload.get("event_type", "se"))
        genome = payload.get("genome", "hg19")
        if event_type not in EVENT_SPECS:
            return jsonify({"success": False, "error": f"Unsupported event type: {event_type}"}), 400

        job_dir = _create_job_dir(Path(app.config["RESULTS_FOLDER"]), prefix="quicktest_")
        job_id = job_dir.name.replace("quicktest_", "")

        with app.jobs_lock:
            app.jobs[job_id] = {
                "status": "queued",
                "event_type": event_type,
                "genome": genome,
                "created": datetime.now(),
                "input_files": {
                    "rmats": str(q["rmats"]),
                    "known_motifs": str(q["known_motifs"]),
                    "motifs": str(q["motifs"]),
                },
                "params": {
                    "window_size": 50,
                    "step_size": 1,
                    "intron_len": 250,
                    "exon_len": 50,
                    "sig_fdr": 0.05,
                    "sig_delta_psi": 0.05,
                    "rbp_label": "QuickTest",
                },
                "output_dir": str(job_dir),
            }

        LOGGER.info("Submitted quick test job %s", job_id)
        app.executor.submit(run_job, job_id)
        return jsonify({"success": True, "job_id": job_id, "output_dir": str(job_dir)}), 200

    @app.post("/api/submit")
    def submit_job():
        try:
            event_type = _normalize_event(request.form.get("event_type", ""))
            genome = request.form.get("genome", "")
            input_type = request.form.get("input_type", "")
            if not event_type or not genome:
                return jsonify({"success": False, "error": "Missing required parameters"}), 400
            if event_type not in EVENT_SPECS:
                return jsonify({"success": False, "error": f"Unsupported event type: {event_type}"}), 400

            job_dir = _create_job_dir(Path(app.config["RESULTS_FOLDER"]))
            job_id = job_dir.name

            input_files: dict[str, str] = {}
            if input_type == "rmats":
                file = request.files.get("rmats_file")
                if not file or not file.filename:
                    return jsonify({"success": False, "error": "Missing rMATS file"}), 400
                path = job_dir / file.filename
                file.save(path)
                input_files["rmats"] = str(path)
            elif input_type == "miso":
                file = request.files.get("miso_file")
                if not file or not file.filename:
                    return jsonify({"success": False, "error": "Missing MISO file"}), 400
                path = job_dir / file.filename
                file.save(path)
                input_files["miso"] = str(path)
            elif input_type == "coordinates":
                for field, key in (("up_file", "up"), ("down_file", "down"), ("bg_file", "bg")):
                    file = request.files.get(field)
                    if not file or not file.filename:
                        return jsonify({"success": False, "error": f"Missing {field}"}), 400
                    path = job_dir / file.filename
                    file.save(path)
                    input_files[key] = str(path)
            else:
                return jsonify({"success": False, "error": "Invalid input type"}), 400

            known = request.files.get("known_motifs_file")
            if not known or not known.filename:
                return jsonify({"success": False, "error": "Missing known motifs file"}), 400
            known_path = job_dir / known.filename
            known.save(known_path)
            input_files["known_motifs"] = str(known_path)

            custom = request.files.get("motifs_file")
            if custom and custom.filename:
                custom_path = job_dir / custom.filename
                custom.save(custom_path)
                input_files["motifs"] = str(custom_path)

            params = {
                "window_size": int(request.form.get("window", 50)),
                "step_size": int(request.form.get("step", 1)),
                "intron_len": int(request.form.get("intron", 250)),
                "exon_len": int(request.form.get("exon", 50)),
                "sig_fdr": float(request.form.get("fdr", 0.05)),
                "sig_delta_psi": float(request.form.get("delta_psi", 0.05)),
                "rbp_label": request.form.get("label", "RBP"),
            }

            with app.jobs_lock:
                app.jobs[job_id] = {
                    "status": "queued",
                    "event_type": event_type,
                    "genome": genome,
                    "created": datetime.now(),
                    "input_files": input_files,
                    "params": params,
                    "output_dir": str(job_dir),
                }

            LOGGER.info("Submitted web job %s", job_id)
            app.executor.submit(run_job, job_id)
            return jsonify({"success": True, "job_id": job_id, "output_dir": str(job_dir)}), 200
        except Exception as exc:
            LOGGER.exception("Submit failed")
            return jsonify({"success": False, "error": str(exc)}), 500

    @app.get("/api/status/<job_id>")
    def get_status(job_id: str):
        job = resolve_job(job_id)
        if not job:
            return jsonify({"success": False, "error": "Job not found"}), 404
        return (
            jsonify(
                {
                    "success": True,
                    "status": job["status"],
                    "event_type": job["event_type"],
                    "genome": job["genome"],
                    "created": job["created"].isoformat(),
                    "error": job.get("error"),
                    "output_dir": job.get("output_dir"),
                    "message": job.get("error") if job.get("status") == "failed" else "",
                }
            ),
            200,
        )

    @app.get("/api/logs/<job_id>")
    def get_logs(job_id: str):
        job = resolve_job(job_id)
        if not job:
            return jsonify({"success": False, "error": "Job not found"}), 404

        output_dir = Path(job["output_dir"])
        log_file = output_dir / "log.motifMap.txt"
        max_lines = request.args.get("tail", default=120, type=int) or 120
        max_lines = min(max(20, max_lines), 400)

        lines: list[str] = []
        if log_file.exists():
            lines = _important_log_lines(_tail_lines(log_file, max_lines * 3))
            lines = lines[-max_lines:]

        return (
            jsonify(
                {
                    "success": True,
                    "job_id": job_id,
                    "status": job["status"],
                    "log_file_exists": log_file.exists(),
                    "lines": lines,
                }
            ),
            200,
        )

    @app.get("/api/results/<job_id>")
    def get_results(job_id: str):
        job = resolve_job(job_id)
        if not job:
            return jsonify({"success": False, "error": "Job not found"}), 404
        if job["status"] not in ("completed", "unknown", "failed"):
            return jsonify({"success": False, "error": "Job not completed"}), 400

        output_dir = Path(job["output_dir"])
        results_file = output_dir / "motif_enrichment_results.txt"

        rows: list[dict[str, Any]] = []
        if results_file.exists():
            with open(results_file, "r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if parts[0] == "RBP" or len(parts) < 5:
                        continue
                    try:
                        rows.append(
                            {
                                "rbp": parts[0],
                                "pval_up_vs_bg": float(parts[1]),
                                "pval_dn_vs_bg": float(parts[2]),
                                "log_pval_up": float(parts[3]),
                                "log_pval_dn": float(parts[4]),
                            }
                        )
                    except Exception:
                        continue
        else:
            # Legacy fallback: derive summary from pVal.*.RNAmap.txt files.
            up_file = output_dir / "pVal.up.vs.bg.RNAmap.txt"
            dn_file = output_dir / "pVal.dn.vs.bg.RNAmap.txt"
            if not up_file.exists() and not dn_file.exists():
                return jsonify({"success": False, "error": "Results file not found"}), 404

            agg: dict[str, dict[str, float]] = {}

            def parse_rnamap(path: Path, key: str) -> None:
                if not path.exists():
                    return
                with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                    _ = handle.readline()
                    for line in handle:
                        line = line.strip()
                        if not line:
                            continue
                        parts = line.split("\t")
                        if len(parts) < 2:
                            continue
                        motif_name = Path(parts[0]).name
                        pvals: list[float] = []
                        for raw in parts[1:]:
                            try:
                                p = float(raw)
                                if p > 0.0:
                                    pvals.append(p)
                            except Exception:
                                continue
                        if not pvals:
                            continue
                        best = min(pvals)
                        row = agg.setdefault(
                            motif_name,
                            {
                                "rbp": motif_name,
                                "pval_up_vs_bg": 1.0,
                                "pval_dn_vs_bg": 1.0,
                                "log_pval_up": 0.0,
                                "log_pval_dn": 0.0,
                            },
                        )
                        row[key] = best

            parse_rnamap(up_file, "pval_up_vs_bg")
            parse_rnamap(dn_file, "pval_dn_vs_bg")

            for row in agg.values():
                row["log_pval_up"] = -math.log10(max(row["pval_up_vs_bg"], 1e-300))
                row["log_pval_dn"] = -math.log10(max(row["pval_dn_vs_bg"], 1e-300))
                rows.append(row)

        rows.sort(key=lambda row: row["pval_up_vs_bg"])
        return jsonify({"success": True, "results": rows[:20]}), 200

    @app.get("/api/jobs")
    def list_jobs():
        limit = request.args.get("limit", default=25, type=int) or 25
        limit = min(max(1, limit), 200)
        results_root = Path(app.config["RESULTS_FOLDER"])
        jobs = []
        if results_root.exists():
            dirs = sorted(
                [d for d in results_root.iterdir() if d.is_dir()],
                key=lambda d: d.stat().st_mtime,
                reverse=True,
            )[:limit]
            for d in dirs:
                raw = d.name
                job_id = raw.replace("quicktest_", "") if raw.startswith("quicktest_") else raw
                jobs.append(
                    {
                        "job_id": job_id,
                        "folder": str(d),
                        "name": raw,
                        "mtime": datetime.fromtimestamp(d.stat().st_mtime).isoformat(),
                        "has_log": (d / "log.motifMap.txt").exists(),
                    }
                )
        return jsonify({"success": True, "jobs": jobs}), 200


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True, host="127.0.0.1", port=5000)
