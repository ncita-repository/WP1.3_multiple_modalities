# Set base image (host OS) (note: Code developed/tested with Python 3.7.10)
FROM python:3.7-slim-buster
#FROM python:3.7-buster # not tried yet

# Set the working directory in the container
WORKDIR C:/Code/WP1.3_multiple_modalities/src/app

## Copy the content of the local src directory to the working directory
#COPY . .

# Copy the directories containing the required modules
COPY . /configs
COPY . /conversion_tools
COPY . /dicom_tools
COPY . /dro_tools
COPY . /fiducials
COPY . /general_tools
COPY . /image_tools
COPY . /io_tools
COPY . /plotting_tools
COPY . /testing
COPY . /xnat_tools

# Copy the main app module
COPY . app.py

# Copy the dependencies file to the working directory
COPY requirements.txt ./

# Install dependencies
#RUN pip install --no-cache-dir -r requirements.txt # The command '/bin/sh -c pip install --no-cache-dir -r requirements.txt' returned a non-zero code: 2
#RUN pip install -r requirements.txt # times out when pip installing SimpleITK
RUN pip install --default-timeout=900 -r requirements.txt

# Command to run on container start
CMD [ "python", "./app.py C:\Code\WP1.3_multiple_modalities\OOP\configs NCITA_TEST_RR2" ]