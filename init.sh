#!/bin/bash

cd backend
sudo node src/index.js postgres root localhost 5432 biomolecules_db &  # start the backend
cd ..
cd emsurfer
ng serve --host 0.0.0.0 --disable-host-check & # start the frontend
cd ..
