'''
 file constant.py
 
 @author Danny Xie 
 
 This file will keep all the static 
 values of the files.
'''
import enum 

class App(enum.Enum):
    '''
    The App enum deal with the app.

    Enum values:
        updater - Updater
    '''
    updater = "Updater" # Represent the application updater.

class TypeErrorEmd(enum.Enum):
    '''
    The TypeErrorEmd enum deal with the type of errors of emd in the updater.

    Enum values:
        DFE - DFE, not a error, default value.
        EDG - Error en la generación del descriptor EMD id.
        ENA - An error occurred (Computed a NaN!).
        ECS - Contour level or std not exist.
        EIN - Error in the insertion of emd.
        EID - Error in the descriptor generation.
        ETS - Error in the time stamp generation.
        EGS - Error in the generation of segments.
        ECE - Error in the calculation of the volume of the EMD.
    '''
    DFE = "DFE" 
    EDG = "Error in the images generation of EMD {0}" 
    ena = "ENA" 
    ECS = "Error in the descriptor generation of EMD {0}.Countour level or std not exist." 
    EIN = "Error in the insertion of EMD {0}" 
    EID = "Error in the descriptor generation of EMD {0}" 
    ETS = "Error in the time stamp generation of EMD {0}" 
    EGS = "Error in the generation of segments EMD {0}" 
    ECE = "Error in the calculation of the volume of the EMD {0}"

class TypeErrorPdb(enum.Enum):
    '''
    The TypeErrorPdb enum deal with the type of errors of pdb in the updater.

    Enum values:
        EEP - Error in execution in PDB.
        ECC - Error in cath complex.
        ECH - Error in the execution cath chain.
        ECD - Error in exeution cath domain.
        EAS - Error in execution atomic structure.
    '''
    EEP = "Error in execution {0} with PDB {1}" 
    ECC = "Error in execution cath complex {0} with complex {1}" 
    ECH = "Error in execution cath chain {0} with chain {1}" 
    ECD = "Error in execution cath domain {0} with domain {1}"
    EAS = "Error in execution atomic structure x emd_entry {0} with atomic_structure {1}, emd_entry {2}" 

class TypeMessage(enum.Enum):
    '''
    The TypeMessage enum deal with the type of messages in the updater.

    Enum values:
        MS1 - Message of the execution start.
        MET - Message of execution time.
        MS2 - Message of execution finished.
    '''
    MS1 = "------------------------------The execution start ------------------------------"
    MET = "Execution time: {0}"
    MS2 = "---------------------------- The execution finished -----------------------------"

class Constant:
    '''
    The Constant class contained some constants used in the updater.

    Constant values:
        formatLog - Format log.
        config    - Configuration format.
    '''
    formatLog = '[%(asctime)s] %(app)s %(file)s [%(type)s] %(message)s  %(error)s'
    config = {'app': App.updater.value, 'file':__file__, 'type': TypeErrorEmd.DFE.value, 'error':''}

class Algorithm_segmentation(enum.Enum):
    '''
    The Algorithm_segmentation enum deal with the algorithm to generate the segments.

    Enum values:
        DEFAULT - default
    '''
    DEFAULT = 'default'