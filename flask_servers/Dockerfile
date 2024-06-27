# Use the official Anaconda image as the base image
FROM python:3.11

# Set the working directory
WORKDIR /app

# Copy the requirements.txt to the container
COPY requirements.txt .

# Install the dependencies
RUN pip install -r requirements.txt
#RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code to the working directory
COPY . .

# Print directory structure
COPY server/db_setup/Homo_sapiens.GRCh37.87.gtf.gz .
COPY server/db_setup/Homo_sapiens.GRCh37.87.gtf.gz /app

# Expose port 5000 to the outside world
EXPOSE 5001

# Define environment variable
ENV FLASK_APP=main.py
ENV FLASK_RUN_HOST=0.0.0.0

# Command to run the application
CMD ["gunicorn", "--bind", "0.0.0.0:5001", "main:app"]