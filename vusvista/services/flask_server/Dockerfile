# Use the official Anaconda image as the base image
FROM python:3.11

# Set the working directory
WORKDIR /flask_server

# add and install requirements
RUN pip install --upgrade pip
COPY requirements.txt .
RUN pip install -r requirements.txt

# add entrypoint.sh
# COPY ./entrypoint.sh .
# RUN chmod +x entrypoint.sh

# Copy the rest of the application code to the working directory
COPY . .

ARG GENE_ANNOTATION_FILE_PATH
COPY ${GENE_ANNOTATION_FILE_PATH} /flask_server

# Expose port 5001 to the outside world
EXPOSE 5001

# Define environment variable
# ENV FLASK_APP=main.py
# ENV FLASK_RUN_HOST=0.0.0.0

# Command to run the application - 5 min timeout
CMD ["gunicorn", "-w", "1", "--bind", "0.0.0.0:5001", "--timeout", "300", "main:app"]
# ENTRYPOINT ["/flask_server/entrypoint.sh"]
#  CMD ["python", "main.py"]
