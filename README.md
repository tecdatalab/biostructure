## Node JS installation

#### Step 1 - Install cURL

```sh
sudo apt install curl
```

#### Step 2 - Add Node JS PPA

```sh
curl -sL https://deb.nodesource.com/setup_10.x | sudo bash -
```

#### Step 3 - Install Node JS

```sh
sudo apt install nodejs
```

#### Step 4 - Install npm

```sh
sudo apt-get install -y npm
```

## Express API installation

#### Step 1 - Change directory to backend

```shell
cd backend
```

#### Step 2 - Install node packages

```shell
npm install
```

#### Step 3 - Run API

```shell
node src/index.js [database username] [database password] [host] [port] [database name]
```

## Angular installation

#### Step 1 - Change directory to emsurfer

```sh
cd emsurfer
```

#### Step 2 - Install node packages

```sh
npm install
```

#### Step 3 - Link angular/cli to npm

```sh
sudo npm link @angular/cli
```

#### Step 4 - Run Angular

```sh
ng serve
```

## Updater installation

#### Step 1 - Update packages

```sh
sudo apt-get update && sudo apt-get upgrade
```

#### Step 2 - Install Python 3.6

```sh
sudo apt-get install python3.6
```

#### Step 3 - Install Pip3

```sh
echo Y | sudo apt-get install python3-pip
```

#### Step 4 - Install Tkinter

```sh
echo Y | sudo apt-get install python3-tk 
```

#### Step 5 - Install OpenGL

```sh
echo Y | sudo apt-get install libglfw3-dev
```

#### Step 6 - Install Libpq-dev python-dev

```sh
sudo apt-get install libpq-dev python-dev
```

#### Step 7 - Install Psycopg2

```sh
sudo pip3 install -U psycopg2
```

#### Step 8 - Install Numpy

```sh
sudo pip3 install -U numpy 
```

#### Step 9 - Install Biopandas

```sh
sudo pip3 install -U biopandas 
```

#### Step 10 - Install Cython

```sh
sudo pip3 install -U Cython 
```

#### Step 11 - Install Pandas

```sh
sudo pip3 install -U pandas
```

#### Step 12 - Install Sklearn

```sh
sudo pip3 install -U sklearn
```

#### Step 13 - Install Triangle

```sh
sudo pip3 install -U triangle
```

#### Step 14 - Install Scipy

```sh
sudo pip3 install -U scipy
```

#### Step 15 - Install Matplotlib

```sh
sudo pip3 install -U matplotlib
```

#### Step 16 - Install Scikit-image

```sh
sudo pip3 install -U scikit-image
```

#### Step 17 - Install Glumpy

```sh
sudo pip3 install -U glumpy
```

#### Step 18 - Install Pyopengl

```sh
sudo pip3 install -U pyopengl
```

#### Step 19 - Install Packaging 

```sh
sudo pip3 install -U packaging 
```

#### Step 20 - Install Appdirs

```sh
sudo pip3 install -U appdirs  
```

#### Step 21 - Install Glfw

```sh
sudo pip3 install -U glfw 
```

#### Step 22 - Install Pyqt5

```sh
sudo pip3 install -U pyqt5
```

#### Step 23 - Install Psycopg2-binary

```sh
sudo pip3 install -U psycopg2-binary
```

#### Step 24 - Install Requests

```sh
sudo pip3 install -U requests
```

#### Step 25 - Permissions of execution to em_volume

```sh
chmod +x db/Updater/generators/em_volume
```

#### Step 26 - Permissions of execution to map2zernike

```sh
chmod +x db/Updater/generators/map2zernike
```

#### Step 27 - Change directory to Main file

```shell
cd biostructure/db/Updater/updater/
```

#### Step 28 - Run Updater

```shell
python3.6 Main.py [log] [initialEMD] [mode] [image] [descriptor]
```
