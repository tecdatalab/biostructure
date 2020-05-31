'''
file constant.py
@author Danny Xie 
This file will keep all the static values of the files.
'''

import enum 

# class TypeErrorEmd enum
# Description: The following class contain all the app that will be using this file.
# By: Danny Xie Li
# Created: 16/05/2020
# Last modification: 16/05/2020

class App(enum.Enum):
    updater = "Updater" # Represent the application updater.

# class TypeErrorEmd enum
# Description: The following class contain all the message of error of the main_emd.py for the logging file.
# By: Danny Xie Li
# Created: 16/05/2020
# Last modification: 16/05/2020

class TypeErrorEmd(enum.Enum):
    DFE = "DFE" # Not a error, default value.
    EDG = "Error in the images generation of EMD {0}" # Error en la generación del descriptor EMD id.
    ena = "ENA" # An error occurred (Computed a NaN!).
    ECS = "Error in the descriptor generation of EMD {0}.Countour level or std not exist." # Contour level or std not exist.
    EIN = "Error in the insertion of EMD {0}" # Error in the insertion.
    EID = "Error in the descriptor generation of EMD {0}" # Error in the descriptor generation.
    ETS = "Error in the time stamp generation of EMD {0}" # Error in the time stamp generation.

# class TypeErrorPdb enum
# Description: The following class contain all the message of error of the main_pdb.py for the logging file.
# By: Danny Xie Li
# Created: 16/05/2020
# Last modification: 16/05/2020

class TypeErrorPdb(enum.Enum):
    EEP = "Error in execution {0} with PDB {1}" # Error in execution in PDB.
    ECC = "Error in execution cath complex {0} with complex {1}" # Error in cath complex.
    ECH = "Error in execution cath chain {0} with chain {1}" # Error in the execution cath chain.
    ECD = "Error in execution cath domain {0} with domain {1}" # Error in exeution cath domain.
    EAS = "Error in execution atomic structure x emd_entry {0} with atomic_structure {1}, emd_entry {2}" # Error in execution atomic structure.

# class TypeMessage enum
# Description: The following class contain all the message of the logging file.
# By: Danny Xie Li
# Created: 16/05/2020
# Last modification: 16/05/2020

class TypeMessage(enum.Enum):
    MS1 = "------------------------------The execution start ------------------------------"
    MET = "Execution time: {0}"
    MS2 = "---------------------------- The execution finished -----------------------------"

# class Constant
# Description: This class contain of the constants of the file
# By: Danny Xie Li
# Created: 16/05/2020
# Last modification: 16/05/2020

class Constant:
    formatLog = '[%(asctime)s] %(app)s %(file)s [%(type)s] %(message)s  %(error)s'
    config = {'app': App.updater.value, 'file':__file__, 'type': TypeErrorEmd.DFE.value, 'error':''}

