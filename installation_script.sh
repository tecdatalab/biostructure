#!/bin/bash

# -> libraries installation

# CURL
sudo apt install curl -y
curl -sL https://deb.nodesource.com/setup_10.x | sudo bash -

# NODE-JS
sudo apt install nodejs
sudo apt-get install -y npm

# ->BACK-END INSTALLATION

# GO INSIDE THE BACK-END FOLDER
cd backend

# INSTALL NODE MODULES
sudo chown -R `whoami` ~/.npm
sudo chown -R `whoami` /usr/local/lib/node_modules
sudo npm install

# GO TO PREVIOUS FOLDER
cd ..

# -> FRONT-END INSTALLATION

# GO TO FRONT-END FOLDER
cd emsurfer

# INSTALL NODE MODULES
sudo chown -R `whoami` ~/.npm
sudo chown -R `whoami` /usr/local/lib/node_modules
npm install
npm audit fix
# ng update @angular/cli @angular/core --force

# LINK ANGULAR CLI TO THE FRONT-END PROJECT
sudo npm link @angular/cli

# GO TO PREVIOUS FOLDER
cd ..

# -> PYTHON LIBRARIES INSTALLATION

# UPDATE & UPGRADE
sudo apt-get update && sudo apt-get upgrade

# PYTHON 3.6
sudo apt-get install python3.6
#sudo apt-get install python-dev
#sudo apt-get install python3-dev

# PIP INSTALL
echo Y | sudo apt-get install python3-pip

# TKINTER INSTALL
echo Y | sudo apt-get install python3-tk

# OpenGL INSTALL
echo Y | sudo apt-get install libglfw3-dev

# Libpq-dev python-dev INSTALL
sudo apt-get install libpq-dev python-dev

#sudo python3 -m pip install --upgrade wheel
#sudo python3 -m pip install --upgrade setuptools

# Psycopg2
sudo python3 -m pip install -U psycopg2

# Numpy
sudo python3 -m pip install -U numpy

# Biopandas
sudo python3 -m pip install -U biopandas

# Cython
sudo python3 -m pip install -U Cython

#Pandas
sudo python3 -m pip install -U pandas

# Sklearn
sudo python3 -m pip install -U sklearn

# Triangle
sudo python3 -m pip install -U triangle

# Scipy
sudo python3 -m pip install -U scipy

# Matplotlib
sudo python3 -m pip install -U matplotlib

# Scikit-image
sudo python3 -m pip install -U scikit-image

# Glumpy
sudo python3 -m pip install -U glumpy

# Pyopengl
sudo python3 -m pip install -U pyopengl

# Packaging
sudo python3 -m pip install -U packaging

# Appdirs
sudo python3 -m pip install -U appdirs

# Glfw
sudo python3 -m pip install -U glfw

# Pyqt5
sudo python3 -m pip install -U pyqt5
#sudo apt-get install python3-pyqt5

# Psycopg2-binary
sudo python3 -m pip install -U psycopg2-binary

# Requests
sudo python3 -m pip install -U requests

# Tabulate
sudo python3 -m pip install tabulate

# Permissions of execution to em_volume
chmod +x db/Updater/generators/em_volume

# Permissions of execution to map2zernike
chmod +x db/Updater/generators/map2zernike

# -> DATA BASE CREATION

# INSTALL POSTGRES (https://www.digitalocean.com/community/tutorials/how-to-install-and-use-postgresql-on-ubuntu-18-04)
sudo apt update
sudo apt install postgresql postgresql-contrib

# RESTORE DATA BASE
psql -h localhost -p 5432 -U postgres < ./db/dev-backup/tecdatalab_DB_backup_esteasv31_15-01-2020.sql

# -> UPDATER RUN

# Change directory to Main file
cd db/Updater/updater/

# Run Updater
#python3.6 Main.py 
