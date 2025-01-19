# Use the official Miniconda image from the Docker Hub
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the environment.yml file
COPY environment.yml /app/

# Create the Conda environment
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "fenicsproject", "/bin/bash", "-c"]

# Install dependencies
COPY . /app/

# Run the Django development server
CMD ["conda", "run", "-n", "fenicsproject", "python", "manage.py", "runserver", "0.0.0.0:8000"]