FROM python:3.10-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Node.js
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/*

# Copy project files
COPY . .

# Install Python dependencies
RUN pip install --no-cache-dir -r web/backend/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Build frontend
WORKDIR /app/web/frontend
RUN npm install
RUN chmod -R 755 node_modules/.bin && npm run build

# Back to root
WORKDIR /app

# Expose port
EXPOSE 8000

# Run backend
CMD ["python", "-m", "uvicorn", "web.backend.app:app", "--host", "0.0.0.0", "--port", "8000"]
