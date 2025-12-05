# Web Interface Deployment Guide

This directory contains a FastAPI backend and React frontend for the RNA Score tool.

## Project Structure

```
web/
├── backend/
│   ├── app.py              # FastAPI application
│   └── requirements.txt     # Python dependencies
└── frontend/
    ├── src/
    │   ├── App.jsx         # Main React component
    │   └── App.css         # Styling
    └── package.json        # Node dependencies
```

## Local Development

### Prerequisites
- Python 3.10+
- Node.js 18+

### Backend Setup

```bash
# Install Python dependencies
pip install -r web/backend/requirements.txt

# Run backend
python -m uvicorn web.backend.app:app --host 0.0.0.0 --port 8000 --reload
```

Backend will be available at `http://localhost:8000`
API docs at `http://localhost:8000/docs`

### Frontend Setup

```bash
# Install Node dependencies
cd web/frontend
npm install

# Run development server
npm run dev
```

Frontend will be available at `http://localhost:5173`

### Using Docker Compose (Recommended)

```bash
docker-compose up
```

This will start both backend and frontend.

## Deployment on Render

### Steps:

1. **Push to GitHub**
   - Create a GitHub repository with your project
   - Push all files including the `Dockerfile`

2. **Create Render Service**
   - Go to https://render.com
   - Click "New +" → "Web Service"
   - Connect your GitHub repository
   - Select branch (e.g., main)

3. **Configure**
   - Name: `rna-score-api`
   - Environment: `Docker`
   - Build Command: (leave default)
   - Start Command: (leave default - uses Dockerfile)
   - Plan: Free tier is fine for testing

4. **Environment Variables** (Optional)
   - Add any needed environment variables in Render dashboard

5. **Deploy**
   - Click "Create Web Service"
   - Render will build and deploy automatically
   - Service will be available at `https://rna-score-api.onrender.com`

### Frontend Deployment (Optional - if separate)

If you prefer to host frontend separately:

1. **Build frontend**
   ```bash
   cd web/frontend
   npm run build
   ```

2. **Deploy to Netlify**
   - Connect GitHub to Netlify
   - Select frontend directory: `web/frontend`
   - Build command: `npm run build`
   - Publish directory: `dist`

3. **Configure API URL**
   - Set environment variable `VITE_API_URL` to your Render backend URL

## API Endpoints

- `POST /api/extract` - Extract distances from structure
- `POST /api/train` - Train scoring function
- `POST /api/score` - Score structures
- `GET /api/job/{job_id}/status` - Check job status
- `GET /api/job/{job_id}/results` - List results
- `GET /api/job/{job_id}/download/{filename}` - Download result file

## Features

- **Extract Tab**: Upload structure files and extract distances
- **Train Tab**: Upload histogram/KDE files to train a scoring model
- **Score Tab**: Score structures using trained models
- **Status Tab**: Monitor job progress and download results

## Notes

- Jobs run in background with unique IDs
- Results stored in `./jobs/` directory (temporary)
- For production, use persistent storage (S3, etc.)
- Large files may need timeout adjustments in Render settings

## Troubleshooting

**Port already in use:**
```bash
lsof -ti:8000 | xargs kill -9
```

**Module not found:**
```bash
pip install -r requirements.txt -r web/backend/requirements.txt
```

**Frontend build issues:**
```bash
cd web/frontend && rm -rf node_modules && npm install
```
