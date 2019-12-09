import os
import http.client as hc
import sys
import json
import time
from tabulate import tabulate
sys.path.append('../')

api_url = "localhost:3000"

def makeDataList(data):
    dataList = []
    counter = 1
    for item in data: 
        temp = []
        temp.append(counter)
        temp.append(item['id_code'])
        temp.append(item['len'])
        temp.append(str(item['class'])+"."+str(item['architecture'])+"."+str(item['topology'])+"."+str(item['homologous']))
        temp.append(item['euclidean_distance'])
        dataList.append(temp)
        counter+=1
    return dataList

def showSimilarMoleculesData(data):
    if data == []:
        print("There's not similar molecules...")
    else:
        print("Information of the similar Molecules")
        print(tabulate(makeDataList(data),
        headers=['Top', 'ID PDB', 'Length', 'CATH', 'Euclidean Distance'],tablefmt="psql",colalign=("center","center","center","center","center",)))
    print("\n")
    print("##############################################################################################################################################")

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
        

def showMoleculeData(dataMolecule, dataCathMolecule):
    print("----------------------------------------------------------------------------------------------------------------------------------------------")
    print("##############################################################################################################################################")
    print("\n")
    print("Information of the molecule")
    print(tabulate([[dataMolecule["id_code"], dataCathMolecule["domain_length"], str(dataCathMolecule["class_number"])+"."+
    str(dataCathMolecule["architecture_number"])+"."+str(dataCathMolecule["topology_number"])+"."+
    str(dataCathMolecule["homologous_superfamily_number"])]],
    headers=['ID PDB', 'Length', 'CATH'],tablefmt="psql",colalign=("center","center","center",)))
    print("\n")

def showErrorMessage():
    print("----------------------------------------------------------------------------------------------------------------------------------------------")
    print("##############################################################################################################################################")
    print("\n")
    print("Sorry but the id of the molecule that you queried, it isn't found in the database, please try with another one...")
    print("\n")
    print("##############################################################################################################################################")

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

def makeQuery(pdbId, optRepre, optTempDB, optCATH, optApplyLenFilter, top):
    dataMolecule = getMolecule(pdbId)
    if not 'message' in dataMolecule:
        dataCathMolecule = getCathMolecule(dataMolecule["id"])
        showMoleculeData(dataMolecule, dataCathMolecule)
        dataSimilarMolecules = getSimilarMolecules(int(dataMolecule["id"]), optRepre-1, optTempDB-1, optCATH-1, optApplyLenFilter-1, top)
        showSimilarMoleculesData(dataSimilarMolecules)
    else:
        showErrorMessage()
    pause()
  
def pause():
    programPause = input("Please, press INTRO to continue...")

def exit():
    for letter in "Thank you for using the system bye...\n":
        sys.stdout.write(letter)
        sys.stdout.flush()
        time.sleep(0.1)

def menuQuery(value):
    print ("""
>====================================================================================<
->Select the surface representation:
\t[1] All atoms
\t[2] Main chain
>====================================================================================<""")
    optRepre = int(input("\tOption select is: "))
    print ("""
>====================================================================================<
->Select the template database:
\t[1] Chain
\t[2] Complex
\t[3] Domain
\t[4] All of above
>====================================================================================<""")
    optTempDB = int(input("\tOption select is: "))
    print ("""
>====================================================================================<
->Select the CATH similarity:
\t[1] None
\t[2] CATH
\t[3] CAT
\t[4] CA
\t[4] C
>====================================================================================<""")
    optCATH = int(input("\tOption select is: "))
    print ("""
>====================================================================================<
->Select the option to apply the length filter:
\t[1] YES
\t[2] NO
>====================================================================================<""")
    optApplyLenFilter = int(input("\tOption select is: "))
    print ("""

>====================================================================================<""")
    top = int(input("\tIndicate the quantity of similar results that you want: "))
    print(""">====================================================================================<""")
    os.system("reset")
    makeQuery(value, optRepre, optTempDB, optCATH, optApplyLenFilter, top)

def menuInput():
    os.system("reset")
    print ("""
Please indicate the next data to make the query:
>====================================================================================<""")
    value = input("\tEnter the id: ")
    print(""">====================================================================================<""")
    menuQuery(value)

def menu():
    while True:
        os.system("reset")
        print(""">====================================================================================<
>===============  CLIENT CONSOLE PYTHON APPLICATION FOR PDB QUERIES  ================<
>====================================================================================<


Hi! Please indicate the option to access to the system or exit:
>====================================================================================<
->Select the option:
\t[1] Make a query
\t[2] Exit
>====================================================================================<""")
        option = input("\tOption selected is: ")
        if option == "1":
            menuInput()
        elif option == "2":
            exit()
            break
        else:
            print("The selected option is invalid, please write a correct one...")

if __name__== "__main__":

    menu()