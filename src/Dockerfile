#syntax=docker/dockerfile:1

# Set base image (host OS) (note: Base OS has Python 3.8.6; Code developed/tested with Python 3.7.10)
#FROM python:3.7-slim-buster
#FROM python:3.7-buster
FROM python:3.8-slim-buster

# Assign environment variables
#ENV PYTHONPATH "${PYTHONPATH}:/"
ENV workdir ./

# Print some info to the console as a check
RUN echo "workdir = $workdir" # ok

# Set the working directory in the container
WORKDIR $workdir

# Copy the content of the local src directory to the working directory
COPY . .

# Install dependencies
RUN pip3 install --upgrade pip
#RUN pip3 install -r requirements.txt # times out when pip installing SimpleITK
RUN pip3 install --default-timeout=900 -r requirements.txt

# Command to run on container start
#CMD ["app.py"]

ENTRYPOINT ["python3"] # usage: docker run copy-roi:0.0.21 -it app.py

#CMD echo docker ps