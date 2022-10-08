# use python3 interpreter.
FROM ubuntu:bionic
FROM python:3
COPY . /bin

# import needed python source files
ADD bin/compare_netcdf.py /

# Update base container install
RUN apt-get update
RUN apt-get upgrade -y

# install command-line tools
RUN apt-get install -y python3-pip locales
RUN python3 -m pip install --upgrade pip

# Ensure locales configured correctly
RUN locale-gen en_US.UTF-8
ENV LC_ALL='en_US.utf8'

# Set python aliases for python3
RUN echo 'alias python=python3' >> ~/.bashrc
RUN echo 'alias pip=pip3' >> ~/.bashrc

# install python dependencies
RUN pip3 install numpy
RUN pip3 install netCDF4 
RUN pip3 install matplotlib
RUN pip3 install xarray
ENTRYPOINT [ "python3", "compare_netcdf.py" ]
