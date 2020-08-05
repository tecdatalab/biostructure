# Updater
#### tecdatalab biostructure project

Updater module of the BioStructure project. For emd and pdb. Can use as an independent program or run it using crontab.

## Technologies
Updater project is created with:
* Python: 3.6 +

## Setup
To run this project, follow the installation of the readme file [README.md](../README.md):
``` bat
$ cd ../db/Updater/updater
```
## Parameters

### Parameters schedule_emd.py
The following are the necessary parameters to run the script `schedule_emd.py`.

| Parameter        | Required      | Description                                                   | 
| -------------    |:-------------:| :-------------                                               |
| -l               | No            | Log file name.                                                | 
| -at              | No            | Attempt EMD coefficient. The number of tries for each EMD.    |
| -ie              | No            | Initial EMD for execution.                                    | 
| -ho              | No            | Time schedule by hour.              | 
| -a               | No            | Amount of emd to execute.              | 
| -m               | Yes           | Execution mode, where "c" is complete and "u" update.         | 
| -i               | Yes           | Generation of images and gif where "Y" is yes and "N" is no.  | 
| -d               | Yes           | Generation of descriptor where "Y" is yes and "N" is no.      | 

### Parameters schedule_pdb.py
The following are the necessary parameters to run the script `schedule_pdb.py`.
| Parameter        | Required      | Description                                                   | 
| -------------    |:-------------:| :-------------                                                |
| -l               | No            | Log file name.                                                | 
| -at              | No            | Attempt EMD coefficient. The number of tries for each EMD.    |
| -ia              | No            | Initial atomic for execution.                                 | 
| -ho              | No            | Time schedule by hour.                                        | 
| -am              | No            | Amount of pdb to execute.                                     | 
| -cco             | Yes           | Generate cath atomic structures for complex.                  | 
| -cc              | Yes           | Generate cath atomic structures for chain.                    | 
| -cd              | Yes           | Generate cath atomic structures for domain.                   | 
| -a               | Yes           | Generation atomic structures.                                 | 
| -ae              | Yes           | Generation connection with pdb and emd.                       | 


## Ouput files
| File             | Description                                  | 
| -------------    |:-------------                                | 
| updater.log      | Log file                                     | 
| emd.csv          | Emd memory and time consumption csv file     |  

## Code Examples
* To run emd script:
``` bat
$ sudo python3 schedule_emd.py -ho 1 -a 2 -m c -i N -d Y -ie 1 -at 2
```

* To run pdb script: 
``` bat
$ sudo python3 schedule_pdb.py -a Y -cco Y -cc Y -cd Y -a Y -ae Y -am 200 -ia 1 -ho 1 -at 2
```
