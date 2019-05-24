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
node src/index.js [database username] [database password] [host] [port]
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
