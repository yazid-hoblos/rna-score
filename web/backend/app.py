"""
FastAPI backend for RNA Scoring Web Interface
"""
import os
import sys
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, List
from fastapi import FastAPI, UploadFile, File, Form, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
import uuid
import json

app = FastAPI(title="RNA Score API", version="1.0.0")

# Path to built frontend (Vite build output)
FRONTEND_DIST = Path(__file__).resolve().parent.parent / "frontend" / "dist"

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve built frontend if present
if FRONTEND_DIST.exists():
    app.mount(
        "/",
        StaticFiles(directory=str(FRONTEND_DIST), html=True),
        name="frontend",
    )

# Base directory for temporary job data
JOBS_DIR = Path("./jobs")
JOBS_DIR.mkdir(exist_ok=True)

# Project root (adjust based on deployment)
PROJECT_ROOT = Path(__file__).parent.parent.parent
SRC_DIR = PROJECT_ROOT / "src"


class JobStatus:
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


def get_job_dir(job_id: str) -> Path:
    """Get or create job directory."""
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(exist_ok=True)
    return job_dir


def save_job_status(job_id: str, status: str, message: str = "", error: str = ""):
    """Save job status to JSON file."""
    job_dir = get_job_dir(job_id)
    status_file = job_dir / "status.json"
    status_file.write_text(json.dumps({
        "status": status,
        "message": message,
        "error": error
    }))


def get_job_status(job_id: str) -> dict:
    """Get job status from JSON file."""
    job_dir = get_job_dir(job_id)
    status_file = job_dir / "status.json"
    if status_file.exists():
        return json.loads(status_file.read_text())
    return {"status": JobStatus.PENDING, "message": "", "error": ""}


@app.get("/api/health")
async def health():
    """Health check endpoint."""
    return {"status": "ok"}


@app.get("/api/job/{job_id}/status")
async def check_status(job_id: str):
    """Check job status."""
    try:
        status = get_job_status(job_id)
        return status
    except Exception as e:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")


@app.post("/api/pipeline")
async def run_pipeline(
    background_tasks: BackgroundTasks,
    train_file: UploadFile = File(...),
    score_file: UploadFile = File(...),
    atom_mode: str = Form("C3'"),
    dist_mode: str = Form("intra"),
    cutoff: float = Form(20.0),
    seq_sep: int = Form(4),
    bin_width: float = Form(1.0),
    method: str = Form("histogram"),
    max_score: float = Form(10.0),
    pseudocount: float = Form(1e-6),
    include_plot: bool = Form(True),
):
    """Run complete pipeline: extract (train) → train → score (test) → plot."""
    job_id = str(uuid.uuid4())
    job_dir = get_job_dir(job_id)
    
    # Save uploaded training and scoring files
    train_path = job_dir / train_file.filename
    with open(train_path, "wb") as f:
        f.write(await train_file.read())
    
    score_path = job_dir / score_file.filename
    with open(score_path, "wb") as f:
        f.write(await score_file.read())
    
    save_job_status(job_id, JobStatus.PENDING, "Pipeline queued for processing")
    
    background_tasks.add_task(
        run_full_pipeline,
        job_id,
        str(train_path),
        str(score_path),
        atom_mode,
        dist_mode,
        cutoff,
        seq_sep,
        bin_width,
        method,
        max_score,
        pseudocount,
        include_plot
    )
    
    return {"job_id": job_id, "status": "queued"}


def run_full_pipeline(job_id: str, train_file: str, score_file: str, atom_mode: str, dist_mode: str,
                      cutoff: float, seq_sep: int, bin_width: float, method: str,
                      max_score: float, pseudocount: float, include_plot: bool):
    """Background task: Run full pipeline (extract train → train → score test → plot)."""
    try:
        job_dir = get_job_dir(job_id)
        print(f"[{job_id}] Starting pipeline. Project root: {PROJECT_ROOT}, SRC_DIR: {SRC_DIR}")
        print(f"[{job_id}] SRC_DIR exists: {SRC_DIR.exists()}")
        print(f"[{job_id}] extract_distances.py exists: {(SRC_DIR / 'extract_distances.py').exists()}")
        
        # Step 1: Extract from training data
        save_job_status(job_id, JobStatus.RUNNING, "Step 1/4: Extracting distances from training data...")
        extract_dir = job_dir / "extracted"
        extract_dir.mkdir(exist_ok=True)
        
        # Determine if input is a list file or a single structure file
        is_list_file = train_file.endswith('.txt')
        print(f"[{job_id}] Training file: {train_file}, is_list: {is_list_file}")
        
        cmd = [
            sys.executable,
            str(SRC_DIR / "extract_distances.py"),
            "--list" if is_list_file else "--pdb", train_file,
            "--atom-mode", atom_mode,
            "--dist-mode", dist_mode,
            "--cutoff", str(cutoff),
            "--seq-sep", str(seq_sep),
            "--bin-width", str(bin_width),
            "--method", method,
            "--out-dir", str(extract_dir)
        ]
        
        print(f"[{job_id}] Running: {' '.join(cmd)}")
        env = os.environ.copy()
        env['PYTHONPATH'] = str(PROJECT_ROOT)
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(PROJECT_ROOT), env=env, timeout=300)
        print(f"[{job_id}] Extract returncode: {result.returncode}")
        if result.stdout:
            print(f"[{job_id}] Extract stdout: {result.stdout}")
        if result.stderr:
            print(f"[{job_id}] Extract stderr: {result.stderr}")
        
        if result.returncode != 0:
            raise Exception(f"Extraction failed: {result.stderr}")
        
        # Step 2: Train
        save_job_status(job_id, JobStatus.RUNNING, "Step 2/4: Training model...")
        training_dir = job_dir / "training_output"
        training_dir.mkdir(exist_ok=True)
        
        cmd = [
            sys.executable,
            str(SRC_DIR / "train.py"),
            "--input-dir", str(extract_dir),
            "--output-dir", str(training_dir),
            "--max-score", str(max_score),
            "--pseudocount", str(pseudocount),
            "--cutoff", str(cutoff),
            "--bin-width", str(bin_width),
            "--method", method
        ]
        
        env = os.environ.copy()
        env['PYTHONPATH'] = str(PROJECT_ROOT)
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(PROJECT_ROOT), env=env, timeout=300)
        if result.returncode != 0:
            raise Exception(f"Training failed: {result.stderr}")
        
        # Step 3: Score test structures
        save_job_status(job_id, JobStatus.RUNNING, "Step 3/4: Scoring test structures...")
        scores_file = job_dir / "scores.csv"
        
        # Determine if input is a list file or a single structure file
        is_list_file = score_file.endswith('.txt')
        
        cmd = [
            sys.executable,
            str(SRC_DIR / "score_structures.py"),
            "--list" if is_list_file else "--pdb", score_file,
            "--tables", str(training_dir),
            "--cutoff", str(cutoff),
            "--seq-sep", str(seq_sep),
            "--output", str(scores_file)
        ]
        
        env = os.environ.copy()
        env['PYTHONPATH'] = str(PROJECT_ROOT)
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(PROJECT_ROOT), env=env, timeout=300)
        if result.returncode != 0:
            raise Exception(f"Scoring failed: {result.stderr}")
        
        # Step 4: Plot
        if include_plot:
            save_job_status(job_id, JobStatus.RUNNING, "Step 4/4: Generating plots...")
            plots_dir = job_dir / "plots"
            plots_dir.mkdir(exist_ok=True)
            
            cmd = [
                sys.executable,
                str(SRC_DIR / "plot_scores.py"),
                "--input-dir", str(training_dir),
                "--output-dir", str(plots_dir)
            ]
            
            env = os.environ.copy()
            env['PYTHONPATH'] = str(PROJECT_ROOT)
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(PROJECT_ROOT), env=env, timeout=300)
            if result.returncode != 0:
                # Don't fail pipeline if plotting fails
                pass
        
        save_job_status(job_id, JobStatus.COMPLETED, "Pipeline complete!")
    
    except subprocess.TimeoutExpired as e:
        error_msg = f"Step timed out after 5 minutes: {str(e)}"
        print(f"[{job_id}] TIMEOUT: {error_msg}")
        save_job_status(job_id, JobStatus.FAILED, "", error=error_msg)
    except Exception as e:
        error_msg = str(e)
        print(f"[{job_id}] ERROR: {error_msg}")
        import traceback
        print(f"[{job_id}] Traceback: {traceback.format_exc()}")
        save_job_status(job_id, JobStatus.FAILED, "", error=error_msg)


@app.get("/api/job/{job_id}/results")
async def list_results(job_id: str):
    """List all result files for a job."""
    job_dir = get_job_dir(job_id)
    results = {
        "job_id": job_id,
        "files": []
    }
    
    if job_dir.exists():
        # Get all files recursively, excluding status.json
        for file in job_dir.rglob("*"):
            if file.is_file() and file.name != "status.json":
                # Store relative path for URL
                rel_path = file.relative_to(job_dir)
                results["files"].append({
                    "name": str(rel_path),
                    "size": file.stat().st_size,
                    "url": f"/api/job/{job_id}/download/{rel_path}"
                })
    
    return results


@app.get("/api/job/{job_id}/download/{file_path:path}")
async def download_file(job_id: str, file_path: str):
    """Download a result file from a job."""
    job_dir = get_job_dir(job_id)
    
    # Prevent directory traversal attacks
    if ".." in file_path or file_path.startswith("/"):
        raise HTTPException(status_code=400, detail="Invalid file path")
    
    full_path = job_dir / file_path
    
    if not full_path.exists() or not full_path.is_file() or "status.json" in file_path:
        raise HTTPException(status_code=404, detail="File not found")
    
    # Determine media type based on file extension
    media_type = 'application/octet-stream'
    filename = full_path.name
    if filename.endswith('.png'):
        media_type = 'image/png'
    elif filename.endswith('.jpg') or filename.endswith('.jpeg'):
        media_type = 'image/jpeg'
    elif filename.endswith('.csv'):
        media_type = 'text/csv'
    elif filename.endswith('.json'):
        media_type = 'application/json'
    elif filename.endswith('.txt'):
        media_type = 'text/plain'
    
    return FileResponse(
        full_path,
        filename=filename,
        media_type=media_type
    )


@app.post("/api/extract")
async def extract_distances(
    background_tasks: BackgroundTasks,
    file: Optional[UploadFile] = File(None),
    atom_mode: str = Form("C3'"),
    dist_mode: str = Form("intra"),
    cutoff: float = Form(20.0),
    seq_sep: int = Form(4),
    bin_width: float = Form(1.0),
    method: str = Form("histogram"),
):
    """Extract distances from uploaded structure file."""
    if not file:
        raise HTTPException(status_code=400, detail="No file uploaded")
    
    job_id = str(uuid.uuid4())
    job_dir = get_job_dir(job_id)
    
    # Save uploaded file
    file_path = job_dir / file.filename
    with open(file_path, "wb") as f:
        f.write(await file.read())
    
    save_job_status(job_id, JobStatus.PENDING, "Job queued for processing")
    
    # Run extraction in background
    background_tasks.add_task(
        run_extract,
        job_id,
        str(file_path),
        atom_mode,
        dist_mode,
        cutoff,
        seq_sep,
        bin_width,
        method
    )
    
    return {"job_id": job_id, "status": "queued"}


def run_extract(job_id: str, pdb_file: str, atom_mode: str, dist_mode: str, 
                cutoff: float, seq_sep: int, bin_width: float, method: str):
    """Background task: run distance extraction."""
    try:
        job_dir = get_job_dir(job_id)
        save_job_status(job_id, JobStatus.RUNNING, "Running distance extraction...")
        
        cmd = [
            sys.executable,
            str(SRC_DIR / "extract_distances.py"),
            "--pdb", pdb_file,
            "--atom-mode", atom_mode,
            "--dist-mode", dist_mode,
            "--cutoff", str(cutoff),
            "--seq-sep", str(seq_sep),
            "--bin-width", str(bin_width),
            "--method", method,
            "--out-dir", str(job_dir / "extracted")
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(PROJECT_ROOT))
        
        if result.returncode == 0:
            # Zip the output
            output_zip = job_dir / "extracted.zip"
            shutil.make_archive(str(output_zip.with_suffix('')), 'zip', str(job_dir / "extracted"))
            save_job_status(job_id, JobStatus.COMPLETED, "Extraction complete")
        else:
            save_job_status(job_id, JobStatus.FAILED, "", error=result.stderr)
    
    except Exception as e:
        save_job_status(job_id, JobStatus.FAILED, "", error=str(e))


@app.post("/api/train")
async def train(
    background_tasks: BackgroundTasks,
    hist_files: List[UploadFile] = File(...),
    max_score: float = Form(10.0),
    pseudocount: float = Form(1e-6),
    cutoff: float = Form(20.0),
    bin_width: float = Form(1.0),
    method: str = Form("histogram"),
):
    """Train scoring function from histogram/KDE files."""
    job_id = str(uuid.uuid4())
    job_dir = get_job_dir(job_id)
    hist_dir = job_dir / "histograms"
    hist_dir.mkdir(exist_ok=True)
    
    # Save uploaded files
    for file in hist_files:
        file_path = hist_dir / file.filename
        with open(file_path, "wb") as f:
            f.write(await file.read())
    
    save_job_status(job_id, JobStatus.PENDING, "Job queued for processing")
    
    background_tasks.add_task(
        run_train,
        job_id,
        str(hist_dir),
        max_score,
        pseudocount,
        cutoff,
        bin_width,
        method
    )
    
    return {"job_id": job_id, "status": "queued"}


def run_train(job_id: str, hist_dir: str, max_score: float, pseudocount: float,
              cutoff: float, bin_width: float, method: str):
    """Background task: run training."""
    try:
        job_dir = get_job_dir(job_id)
        output_dir = job_dir / "training_output"
        
        save_job_status(job_id, JobStatus.RUNNING, "Training scoring function...")
        
        cmd = [
            sys.executable,
            str(SRC_DIR / "train.py"),
            "--input-dir", hist_dir,
            "--output-dir", str(output_dir),
            "--max-score", str(max_score),
            "--pseudocount", str(pseudocount),
            "--cutoff", str(cutoff),
            "--bin-width", str(bin_width),
            "--method", method
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(PROJECT_ROOT))
        
        if result.returncode == 0:
            output_zip = job_dir / "training_output.zip"
            shutil.make_archive(str(output_zip.with_suffix('')), 'zip', str(output_dir))
            save_job_status(job_id, JobStatus.COMPLETED, "Training complete")
        else:
            save_job_status(job_id, JobStatus.FAILED, "", error=result.stderr)
    
    except Exception as e:
        save_job_status(job_id, JobStatus.FAILED, "", error=str(e))


@app.post("/api/score")
async def score(
    background_tasks: BackgroundTasks,
    files: List[UploadFile] = File(...),
    tables_dir: Optional[UploadFile] = File(None),
    cutoff: float = Form(20.0),
    seq_sep: int = Form(4),
):
    """Score structures."""
    job_id = str(uuid.uuid4())
    job_dir = get_job_dir(job_id)
    structures_dir = job_dir / "structures"
    structures_dir.mkdir(exist_ok=True)
    
    # Save structure files
    for file in files:
        file_path = structures_dir / file.filename
        with open(file_path, "wb") as f:
            f.write(await file.read())
    
    # Use default tables or uploaded
    tables_path = str(SRC_DIR.parent / "training_output")
    
    save_job_status(job_id, JobStatus.PENDING, "Job queued for processing")
    
    background_tasks.add_task(
        run_score,
        job_id,
        str(structures_dir),
        tables_path,
        cutoff,
        seq_sep
    )
    
    return {"job_id": job_id, "status": "queued"}


def run_score(job_id: str, structures_dir: str, tables_dir: str, cutoff: float, seq_sep: int):
    """Background task: run scoring."""
    try:
        job_dir = get_job_dir(job_id)
        output_file = job_dir / "scores.csv"
        
        save_job_status(job_id, JobStatus.RUNNING, "Scoring structures...")
        
        cmd = [
            sys.executable,
            str(SRC_DIR / "score_structures.py"),
            "--folder", structures_dir,
            "--tables", tables_dir,
            "--cutoff", str(cutoff),
            "--seq-sep", str(seq_sep),
            "--output", str(output_file)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(PROJECT_ROOT))
        
        if result.returncode == 0:
            save_job_status(job_id, JobStatus.COMPLETED, "Scoring complete")
        else:
            save_job_status(job_id, JobStatus.FAILED, "", error=result.stderr)
    
    except Exception as e:
        save_job_status(job_id, JobStatus.FAILED, "", error=str(e))


# SPA fallback: serve frontend for non-API routes
@app.get("/{full_path:path}")
async def spa_fallback(full_path: str):
    """Serve index.html for any non-API route (supports client-side routing)."""
    if full_path.startswith("api/"):
        raise HTTPException(status_code=404, detail="Not Found")
    index_file = FRONTEND_DIST / "index.html"
    if index_file.exists():
        return FileResponse(index_file)
    raise HTTPException(status_code=404, detail="Not Found")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
