#!/bin/bash

cd backend
gnome-terminal -- node src/index.js postgres root localhost 5432 biomolecules_db  # start the backend
cd ..
cd emsurfer
gnome-terminal -- ng serve --host 0.0.0.0 --disable-host-check  # start the frontend
cd ..
