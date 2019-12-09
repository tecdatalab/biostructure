# Guide for Python PDB console client

This is a little guide of how you can use the console client to make queries to the PDB structures usign the Python program. The guide consist in 3 parts:
1. How to use
2. Code documentation
3. Extra

## How to use

This section consists to explain how to use the Python program. The first part is an explanation related with some libraries that you have to install to run the program and the second part is a demonstration step by step of how interact with the program.

### Program Installation

Something that is a restriction to use the program is that you need Python3 to compile and run the program, also you need a code editor and install a library. To see more information about some code editors and how to install them, go to the [**Extra**](#extra) section


#### Libraries

Then you have to install the next library to run de Python program:

##### Tabulate

To install the library tabulate, you just need to pen the terminal and run the next command

```
sudo pip3 install tabulate
```
Is necessary to install this library, because is the one the show the information in a format that can be changed by the user.

### Demostration

To run the client, you have to follow the next steps:

- You need to go to the folder where the document is, using the terminal

![](https://github.com/tecdatalab/biostructure/blob/master/client/images/Screenshot%20from%202019-11-25%2022-21-37.png)

- To run the file you need to run the next command in a terminal: ```python3 console_client_PDB.py```, after that the system displays a menu

![](https://github.com/tecdatalab/biostructure/blob/master/client/images/Screenshot%20from%202019-11-25%2022-25-51.png)

- Enter the key 1 to make a query, in case of enter the key 2 the program ends.

![](https://github.com/tecdatalab/biostructure/blob/master/client/images/Screenshot%20from%202019-11-25%2023-06-15.png)

- Enter the id of the molecule that you want to query and example of an id can be *101mA*, and enter the key numbers of the options that you want to apply to the query.

![](https://github.com/tecdatalab/biostructure/blob/master/client/images/Screenshot%20from%202019-11-25%2023-11-00.png)

- Then you can see the results of the query, to continue using the system you need to press the key INTRO and the system display the menu again

![](https://github.com/tecdatalab/biostructure/blob/master/client/images/Screenshot%20from%202019-11-25%2023-31-08.png)
  
The result depends of the data that is in the database, so if you don't have a lot of data, it's possible that your result will be little, but if you have the complete database, you must have the complete result of the query.

## Code documentation

This section consists to explain imports of libraries, variables and functions of the Python program. The comments of the functions are going to be: a general description, the inputs, the outputs and if there are restrictions, the restrictions of the function.

### Imports

The program have some imports of some libraries that are needed to execute itself, those are:

- import os: import the operative system
- import http.client as hc: import for the library http
- import sys: import the system
- import json: import the library Json
- import time: import the time
- from tabulate import tabulate: import the library tabulate

### Variables

The variables of the program are:

- sys.path.append('../'): this is to define the folder where the program is located as the default path
- api_url = "localhost:3000": this is the url where the server is located. If the server is located in other URL, you just need to change the actual URL to the new one. Something to consider when you want to change the URL is that you just copy the URL after the "http://", for example id the URL is "http://localhost:3000", the new URL is "localhost:3000".


### Main functions

#### makeQuery

This is the main function of the program cause is the function that controls and manages the order of how the others functions call. This function receives 6 inputs and doesn't return outputs.

```python
def makeQuery(pdbId, optRepre, optTempDB, optCATH, optApplyLenFilter, top):
    dataMolecule = getMolecule(pdbId)
    dataCathMolecule = getCathMolecule(dataMolecule["id"])
    showMoleculeData(dataMolecule, dataCathMolecule)
    dataSimilarMolecules = getSimilarMolecules(int(dataMolecule["id"]), optRepre-1, optTempDB-1, optCATH-1, optApplyLenFilter-1, top)
    showSimilarMoleculesData(dataSimilarMolecules)
    pause()
```

##### Inputs

The 6 inputs of the funtion are:

- pdbId
  - kind of data: string.
  - value: '100md'
  - meaning: id of the molecule that you want to query in to the database.
  - restrictions: can't be null.
- optRepre
  - kind of data: integer.
  - value: [1, 2]
  - meaning: kind of surface representation to the query.
  - restrictions: can't be null.
- optTempDB
  - kind of data: integer.
  - value: [1 ,2, 3, 4]
  - meaning: template database to the query.
  - restrictions: can't be null.
- optCATH
  - kind of data: integer.
  - value: [1, 2]
  - meaning: applies the cath similarity filter to the query.
  - restrictions: can't be null.
- optApplyLenFilter
  - kind of data: integer.
  - value: [1, 2]
  - meaning: applies the length filter to the query.
  - restrictions: can't be null.
- top
  - kind of data: integer.
  - value: [25, 50, 100, 500, 1000]
  - meaning: quantity of similar molecules that you want to return on the top.
  - restrictions: can't be null.

#### getMolecule

This is the function that obtains the general information of the molecule that is querying to the database. To obtain the information the function uses a HTTP GET Request to the Node Server and adds to the string request the id of the molecule. This functions receives only one input and returns dictionary with the data.

```python
def getMolecule(pdbId):
    try:
        conn = hc.HTTPConnection(api_url)
        stringRequest = "/search/pdb/{}".format(pdbId)
        conn.request("GET", stringRequest)
        request = conn.getresponse()
        result = request.read()
        dictData = json.loads(result.decode('utf-8'))
        conn.close()
        return dictData
    except hc.error as identifier:
        conn.close()
```

##### Inputs

The input of the funtion is:

- pdbId
  - kind of data: string.
  - value: '101mA'
  - meaning: id of the molecule that you want to query in to the database.
  - restrictions: can't be null.

##### Outputs

The output of the funtion is:

- dicData
  - kind of data: dictionary.
  - value: ``` dicData = {'id': 8,'id_code': '101mA', 'numbers_descriptor': [121, 1.24420208483934, 1.60521561000496, ...  ,0.0154185779392719, 0.355090536177158, 1.27513810154051, 0.89794635027647], 'sequence': 'MET\nVAL\nLEU\nSER\nGLU\nGLY\nGLU\nTRP\nGLN\nLEU\nVAL\nLEU', 'png_img_3d': None, 'gif_img_3d': None} ```
  - meaning: general data of the molecule.
  - restrictions: can't be null.  

#### getCathMolecule

This is the function that obtains the CATH information of the molecule that is querying to the database. To obtain the information the function uses a HTTP GET Request to the Node Server and adds to the string request the id key of the molecule on the database. This functions receives only one input and returns dictionary with the data.

```python
def getCathMolecule(pdbId):
    try:
        conn = hc.HTTPConnection(api_url)
        stringRequest = "/search/{}".format(pdbId)
        conn.request("GET", stringRequest)
        request = conn.getresponse()
        result = request.read()
        dictData = json.loads(result.decode('utf-8'))
        conn.close()
        return dictData
    except hc.error as identifier:
        conn.close()
```

##### Inputs

The input of the funtion is:

- pdbId
  - kind of data: integer.
  - value: 1
  - meaning: id key of the molecule on the database that you want to query in to the database.
  - restrictions: can't be null.

##### Outputs

The output of the funtion is:

- dicData
  - kind of data: dictionary.
  - value: ``` dicData = {} ```
  - meaning: CATH data of the molecule.
  - restrictions: can't be null.  

#### getSimilarMolecules

This is the function that obtains the information of the molecules that are similar to the molecule that was querying using different parameters. To obtain the information the function uses a HTTP GET Request to the Node Server and adds the parameters to make the query. This functions receives 6 inputs and returns dictionary with the data.

```python
def getSimilarMolecules(id, optRepre, optTempDB, optCATH, optApplyLenFilter, top):
    try:
        conn = hc.HTTPConnection(api_url)
        stringRequest = "/search/getResultsPDB/{}/{}/{}/{}/{}/{}".format(id, optRepre, optTempDB, optCATH, optApplyLenFilter, top)
        conn.request("GET", stringRequest)
        request = conn.getresponse()
        result = request.read()
        listData = json.loads(result.decode('utf-8'))
        conn.close()
        return listData[0]
    except hc.error as identifier:
        conn.close()
```

##### Inputs

The 6 inputs of the funtion are:

- pdbId
  - kind of data: string.
  - value: '100md'
  - meaning: id of the molecule that you want to query in to the database.
  - restrictions: can't be null.
- optRepre
  - kind of data: integer.
  - value: [1, 2]
  - meaning: kind of surface representation to the query.
  - restrictions: can't be null.
- optTempDB
  - kind of data: integer.
  - value: [1 ,2, 3, 4]
  - meaning: template database to the query.
  - restrictions: can't be null.
- optCATH
  - kind of data: integer.
  - value: [1, 2]
  - meaning: applies the cath similarity filter to the query.
  - restrictions: can't be null.
- optApplyLenFilter
  - kind of data: integer.
  - value: [1, 2]
  - meaning: applies the length filter to the query.
  - restrictions: can't be null.
- top
  - kind of data: integer.
  - value: [25, 50, 100, 500, 1000]
  - meaning: quantity of similar molecules that you want to return on the top.
  - restrictions: can't be null.

##### Outputs

The output of the funtion is:

- listData
  - kind of data: list of dictionaries.
  - value: ``` listData = [{'type_id': 3, 'id': 154, 'id_code': '11asA00', 'id_parent': 151, 'len': 42, 'class': 5, 'architecture': 1, 'topology': 1, 'homologous': 1, 'factor': '3.5000000000000000', 'euclidean_distance': 7.84133824762437},  ... ,{'type_id': 3, 'id': 63, 'id_code': '107mA00', 'id_parent': 61, 'len': 19, 'class': 5, 'architecture': 1, 'topology': 1, 'homologous': 1, 'factor': '1.5833333333333333', 'euclidean_distance': 8.20006168153436}]```
  - meaning: data of the similar molecules.
  - restrictions: can't be null.  

## Extra

### Code editors & Text editors

In my case I use Visual Studio to edit the program but there are some others code or text editors to do that task. Now, I'm going to show how to install some of these:

#### Visual Code

##### Windows

Steps
1. Go to this link: https://code.visualstudio.com/download
2. Click the option that says Windows
3. Execute the downloaded file
4. Open VSC

##### Linux

Steps
1. Open the terminal and run those commands

```
sudo apt update
sudo apt install software-properties-common apt-transport-https wget
wget -q https://packages.microsoft.com/keys/microsoft.asc -O- | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://packages.microsoft.com/repos/vscode stable main"
sudo apt update
sudo apt install code
```
2. Open VSC

#### Sublime Text 3

##### Windows

Steps
1. Go to this link: https://www.sublimetext.com/3
2. Click the option that says Windows
3. Execute the downloaded file
4. Open Sublime Text

##### Linux

Steps
1. Open the terminal and run those commands

```
sudo apt update
sudo add-apt-repository ppa:webupd8team/sublime-text-3
sudo apt-get update
sudo apt-get install sublime-text-installer
```
2. Open Sublime Text

#### Notepad++

##### Windows

Steps
1. Go to this link: https://notepad-plus-plus.org/downloads/
2. Click the last release7.8.1
3. Click the option that says Windows
4. Execute the downloaded file
5. Open Notepad++

##### Linux

Steps
1. Open the terminal and run those commands

```
sudo apt list
sudo apt-get install snapd snapd-xdg-open
sudo snap install notepad-plus-plus
```
2. Open Notepad++

#### PyCharm

##### Windows

Steps
1. Go to this link: https://www.jetbrains.com/pycharm/download/#section=windows
2. Click Community version
3. Execute the downloaded file
4. Open JetBrains

##### Linux

Steps
1. Open the terminal and run those commands

```
sudo add-apt-repository ppa:mystic-mirage/pycharm
sudo apt-get update
sudo apt-get install pycharm-community
```
2. Open JetBrains
