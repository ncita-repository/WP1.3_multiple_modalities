 #syntax=docker/dockerfile:1

# Set base image (host OS) (note: Base OS has Python 3.8.6; Code developed/tested with Python 3.7.10)
#FROM python:3.7-slim-buster
#FROM python:3.7-buster
FROM python:3.8-slim-buster

# Assign environment variables
#ENV PYTHONPATH "${PYTHONPATH}:C:\Code\WP1.3_multiple_modalities\src"
#ENV PYTHONPATH "${PYTHONPATH}:/"
#ENV workdir /src
#ENV workdir ./src
ENV workdir ./
#ENV cfgDir $workdir/configs
#ENV runID # expects a value -- not required for build?

# Print some info to the console as a check
#RUN echo "workdir = ${workdir}" # wrong?
#RUN echo "cfgDir = ${cfgDir}" # wrong?
RUN echo "workdir = $workdir" # ok
#RUN echo "cfgDir = $cfgDir" # ok
#RUN echo "workdir = "$workdir # ok
#RUN echo "cfgDir = "$cfgDir # ok


# Set the working directory in the container
#WORKDIR C:/Code/WP1.3_multiple_modalities/src
#WORKDIR /
WORKDIR $workdir


# Copy the content of the local src directory to the working directory
COPY . .


# Add the app module
#ADD app.py .

# Copy the dependencies file to the working directory
#COPY requirements.txt ./
#COPY requirements.txt requirements.txt
#COPY requirements.txt .

# Install dependencies
RUN pip3 install --upgrade pip
#RUN pip3 install --no-cache-dir -r requirements.txt # The command '/bin/sh -c pip install --no-cache-dir -r requirements.txt' returned a non-zero code: 2
#RUN pip3 install -r requirements.txt # times out when pip installing SimpleITK
RUN pip3 install --default-timeout=900 -r requirements.txt


# Copy the cfgDir (which may or may not be /configs
#COPY . $cfgDir

# Copy the directories containing the required modules
#COPY . /configs
#COPY . /conversion_tools
#COPY . /dicom_tools
#COPY . /dro_tools
#COPY . /fiducials
#COPY . /general_tools
#COPY . /image_tools
#COPY . /io_tools
#COPY . /plotting_tools
#COPY . /sample_DROs
#COPY . /testing
#COPY . /xnat_tools

# Copy the main app module
#COPY . app.py


# Command to run on container start
#CMD [ "python3", "./app.py C:/Code/WP1.3_multiple_modalities/src/configs NCITA_TEST_RR2" ]
#CMD [ "python3", "./app.py ./configs NCITA_TEST_RR2" ]
#CMD [ "python3", "app.py /configs NCITA_TEST_RR2" ]
#CMD [ "python3", "C:/Code/WP1.3_multiple_modalities/src/app.py C:/Code/WP1.3_multiple_modalities/src/configs NCITA_TEST_RR2" ]
#CMD [ "python3", ".app.py ./configs NCITA_TEST_RR2" ]
#CMD [ "python3", ".app.py ./configs NCITA_TEST_RR2" ]
#CMD ["app.py ./configs NCITA_TEST_RR2"]

#CMD ["app.py"]

# Update the configuration files (so that the ENV workdir can be used in directory paths):
#CMD ["config.py ./configs"] # python3: can't open file 'config.py ./configs': [Errno 2] No such file or directory
#CMD ["config.py ${cfgDir}"] # try this?
#CMD ["./config.py ./configs"] 

ENTRYPOINT ["python3"] # usage: docker run copy-roi:0.0.21 -it app.py ./configs NCITA_test_RR2
# FileNotFoundError: [Errno 2] No such file or directory: 'C:\\Code\\WP1.3_multiple_modalities\\src\\inputs\\sample_DROs/sample_spatial_dro.dcm'

#ENTRYPOINT ["python3", ".app.py"] # usage: docker run --rm copy-roi:0.0.22 ./configs NCITA_test_RR2
# python3: can't open file '.app.py': [Errno 2] No such file or directory

#ENTRYPOINT ["python3", "./app.py"] # usage: docker run --rm copy-roi:0.0.22 ./configs NCITA_test_RR2
# FileNotFoundError: [Errno 2] No such file or directory: 'C:\\Code\\WP1.3_multiple_modalities\\src\\inputs\\sample_DROs/sample_spatial_dro.dcm'


#CMD echo docker ps