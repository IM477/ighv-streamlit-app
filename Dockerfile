FROM python:3.9-slim

# Install required tools
RUN apt-get update && \
    apt-get install -y wget unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and install CAP3
RUN wget http://doua.prabi.fr/archive/cap3/cap3.linux.x86_64.tar && \
    tar -xvf cap3.linux.x86_64.tar && \
    mv cap3 /usr/local/bin && \
    chmod +x /usr/local/bin/cap3 && \
    rm cap3.linux.x86_64.tar readme.txt

# Set working directory
WORKDIR /app

# Copy app files
COPY . .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Expose Streamlit port
EXPOSE 8501

# Start Streamlit app
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
