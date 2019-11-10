'''
Created on 30 mar. 2019

@author: luis98
'''
from ftplib import FTP
from mycellanic_classes.pdb_result import Pdb_result
from classes.cath_atomic_structure import Cath_atomic_structure
from mycellanic_classes.cath_domain_boundaries_result import Cath_domain_boundaries_result
from xml.dom import minidom
from classes.atomic_structure_x_emd_entry import Atomic_structure_x_emd_entry
import datetime
import urllib.request
import os
from operator import add


def max_datetime(date1, date2):
    if date1 is None:
        return date2
    elif date2 is None:
        return date1
    elif date1 > date2:
        return date1
    else:
        return date2


class FTP_connection(object):
    '''
    classdocs
    '''
    __ftp = None

    def __init__(self, ftp_server='ftp.ebi.ac.uk', port=None):
        if port is None:
            self.__ftp = FTP()
            self.__ftp.connect(ftp_server)
        else:
            self.__ftp = FTP()
            self.__ftp.connect(ftp_server, port)

    def init_connection(
            self,
            initial_directory='/pub/databases/emdb/structures/'):
        self.__ftp.login()
        if initial_directory is not None:
            # change into "structures" directory
            self.__ftp.cwd(initial_directory)

    def close_connection(self):
        self.__ftp.close()

    def sortKey(self, val):
        return int(val)

    def get_all_emds_id(self, initialEMD, finalEMD):
        files = []
        result = []
        #elements = [172,174,306,307,436,616,1508,1571,1653,1662,1823,1824,1826,1906,2061,2062,2101,2272,2273,2362,2363,2432,2489,2610,2641,2708,3036,3094,3095,3122,3771,3792,3820,3821,3868,4065,4263,4376,4377,4419,4421,4422,4433,4445,4446,4447,4461,4483,4484,4486,4490,4491,4492,4633,5003,5160,5384,5485,5486,5538,5539,5582,5584,5678,5695,5696,5724,5728,5754,5764,5947,5992,6034,6035,6037,6038,6148,6213,6214,6371,6374,6375,6376,6377,6407,6444,6471,6472,6473,6619,6620,6737,6846,6854,6855,6907,6914,6948,6968,6969,6976,7034,7047,7110,7116,7117,7125,7126,7135,7136,7138,7139,7140,7141,7142,7143,7144,7145,7146,7147,7148,7150,7151,7152,7153,7154,7315,7316,7317,7472,7543,7589,7599,7623,7624,7625,7627,7628,7629,7630,8116,8147,8157,8304,8314,8320,8328,8332,8333,8334,8344,8471,8485,8606,8655,8656,8661,8689,8703,8733,8734,8753,8754,8755,8760,8761,8762,8763,8846,8854,9012,9022,9198,9199,9200,9260,9261,9366,9367,9564,9565,9665,9666,9667,9864,9865,9866,9867,9868,9869,9879,9920,9921,20061,20067,20068,20091,20195,20197,20199,20206]
        self.__ftp.retrlines("NLST", files.append)
        for i in range(len(files)):
            # if int(files[i].replace("EMD-", "")) in elements:
            #    continue
            if int(files[i].replace("EMD-", "")) < int(initialEMD):
                continue
            if int(files[i].replace("EMD-", "")) > float(finalEMD):
                break
            try:
                self.__ftp.voidcmd(
                    "MDTM EMD-{0}/map/emd_{0}.map.gz".format(files[i].replace("EMD-", "")))
                result.append(files[i].replace("EMD-", ""))
            except BaseException:
                pass
        '''for i in range(0,4):
            try:
                self.__ftp.voidcmd("MDTM EMD-{0}/map/emd_{0}.map.gz".format(files[i].replace("EMD-", "")))
                result.append(files[i].replace("EMD-", ""))
            except:
                pass'''
        result.sort(key=self.sortKey)
        return result

    def get_emds_higher_than_date(self, date, initialEMD, finalEMD):
        emds_id = self.get_all_emds_id(initialEMD, finalEMD)
        result = []
        for emd in emds_id:
            time_xml = None
            time_image = None
            time_map = None
            try:
                time_xml = self.__ftp.voidcmd(
                    "MDTM EMD-{0}/header/emd-{0}.xml".format(emd))
                time_xml = datetime.datetime.strptime(
                    time_xml[4:], "%Y%m%d%H%M%S")
            except BaseException:
                pass
            try:
                time_image = self.__ftp.voidcmd(
                    "MDTM EMD-{0}/images/emd_{0}.png".format(emd))
                time_image = datetime.datetime.strptime(
                    time_image[4:], "%Y%m%d%H%M%S")
            except BaseException:
                pass
            try:
                time_map = self.__ftp.voidcmd(
                    "MDTM EMD-{0}/map/emd_{0}.map.gz".format(emd))
                time_map = datetime.datetime.strptime(
                    time_map[4:], "%Y%m%d%H%M%S")
            except BaseException:
                pass

            max_value = max_datetime(
                max_datetime(
                    time_xml,
                    time_image),
                time_map)

            if max_value > date:
                result.append(emd)

        return result

    def get_ftp(self):
        return self.__ftp

    def set_ftp(self, value):
        self.__ftp = value

    def del_ftp(self):
        del self.__ftp

    def get_all_pdb(
            self,
            url="http://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt"):
        urllib.request.urlretrieve(url, 'pdb_entry_type.txt')
        file = open("pdb_entry_type.txt", "r")
        result = []
        for line in file:
            fields = line.split("\t")
            data1 = fields[0]
            data2 = fields[1]
            data3 = fields[2]
            result.append(Pdb_result(data1, data2, data3))
        file.close()
        os.remove("pdb_entry_type.txt")
        return result
    
    def get_all_cath_complex(
            self,
            url="http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-chain-list.txt"):
        urllib.request.urlretrieve(url, 'cath-chain-list.txt')
        file = open("cath-chain-list.txt", "r")
        dic_complex = {}
        result = []
        for line in file:
            fields = line.split()
            if fields[0][0] == "#":
                continue
            
            if len(fields)<12:
                data1 = fields[0]
                data2 = fields[1]
                data3 = fields[2]
                data4 = fields[3]
                data5 = fields[4]
                data6 = fields[5]
                data7 = fields[6]
                data8 = fields[7]
                data9 = fields[8][0]
                data10 = fields[8][1:-1]
                data11 = fields[9]
                data12 = fields[10]
            else:
                data1 = fields[0]
                data2 = fields[1]
                data3 = fields[2]
                data4 = fields[3]
                data5 = fields[4]
                data6 = fields[5]
                data7 = fields[6]
                data8 = fields[7]
                data9 = fields[8]
                data10 = fields[9]
                data11 = fields[10]
                data12 = fields[11]
                
            var = [int(data2),int(data3),int(data4),int(data5),int(data6),int(data7),int(data8),int(data9),int(data10),int(data11),float(data12)]
            if dic_complex.get(data1[0:4]) == None:
                dic_complex[data1[0:4]] = var
            else:
                dic_complex[data1[0:4]] =  list( map(add, dic_complex[data1[0:4]], var) )
                

        for key, value in dic_complex.items():
            result.append(
                Cath_atomic_structure(
                    None,
                    key,
                    value[0],
                    value[1],
                    value[2],
                    value[3],
                    value[4],
                    value[5],
                    value[6],
                    value[7],
                    value[8],
                    value[9],
                    value[10]))
            file.close()
        os.remove("cath-chain-list.txt")
        return result

    def get_all_cath_chain(
            self,
            url="http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-chain-list.txt"):
        urllib.request.urlretrieve(url, 'cath-chain-list.txt')
        file = open("cath-chain-list.txt", "r")
        result = []
        for line in file:
            fields = line.split()
            if fields[0][0] == "#":
                continue
            if len(fields)<12:
                data1 = fields[0]
                data2 = fields[1]
                data3 = fields[2]
                data4 = fields[3]
                data5 = fields[4]
                data6 = fields[5]
                data7 = fields[6]
                data8 = fields[7]
                data9 = fields[8][0]
                data10 = fields[8][1:-1]
                data11 = fields[9]
                data12 = fields[10]
            else:
                data1 = fields[0]
                data2 = fields[1]
                data3 = fields[2]
                data4 = fields[3]
                data5 = fields[4]
                data6 = fields[5]
                data7 = fields[6]
                data8 = fields[7]
                data9 = fields[8]
                data10 = fields[9]
                data11 = fields[10]
                data12 = fields[11]
            result.append(
                Cath_atomic_structure(
                    None,
                    data1,
                    data2,
                    data3,
                    data4,
                    data5,
                    data6,
                    data7,
                    data8,
                    data9,
                    data10,
                    data11,
                    data12))
        file.close()
        os.remove("cath-chain-list.txt")
        return result

    def get_all_cath_domain_boundarie(
            self,
            url="http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-boundaries-seqreschopping.txt"):
        urllib.request.urlretrieve(
            url, 'cath-domain-boundaries-seqreschopping.txt')
        file = open("cath-domain-boundaries-seqreschopping.txt", "r")
        result = []
        for line in file:
            fields = line.split()
            data1 = fields[0]
            temp = fields[1].split(",")
            data2 = []
            for i in temp:
                temp_range = i.split("-")
                data2.append([int(temp_range[0]), int(temp_range[1])])
            result.append(Cath_domain_boundaries_result(data1, data2))
        file.close()
        os.remove("cath-domain-boundaries-seqreschopping.txt")
        return result

    def get_all_cath_domain_boundarie_dic(
            self,
            url="http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-boundaries-seqreschopping.txt"):
        urllib.request.urlretrieve(
            url, 'cath-domain-boundaries-seqreschopping.txt')
        file = open("cath-domain-boundaries-seqreschopping.txt", "r")
        dic = {}
        for line in file:
            fields = line.split()
            data1 = fields[0]
            temp = fields[1].split(",")
            data2 = []
            for i in temp:
                temp_range = i.split("-")
                data2.append([int(temp_range[0]), int(temp_range[1])])
            element = Cath_domain_boundaries_result(data1, data2)
            if dic.get(data1[0:4]) is None:
                dic[data1[0:4]] = [element]
            else:
                dic[data1[0:4]].append(element)

        file.close()
        os.remove("cath-domain-boundaries-seqreschopping.txt")
        return dic

    def get_all_cath_domain(
            self,
            url="http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt"):
        urllib.request.urlretrieve(url, 'cath-domain-list.txt')
        file = open("cath-domain-list.txt", "r")
        result = []
        for line in file:
            fields = line.split()
            if fields[0][0] == "#":
                continue
            if len(fields)<12:
                data1 = fields[0]
                data2 = fields[1]
                data3 = fields[2]
                data4 = fields[3]
                data5 = fields[4]
                data6 = fields[5]
                data7 = fields[6]
                data8 = fields[7]
                data9 = fields[8][0]
                data10 = fields[8][1:-1]
                data11 = fields[9]
                data12 = fields[10]
            else:
                data1 = fields[0]
                data2 = fields[1]
                data3 = fields[2]
                data4 = fields[3]
                data5 = fields[4]
                data6 = fields[5]
                data7 = fields[6]
                data8 = fields[7]
                data9 = fields[8]
                data10 = fields[9]
                data11 = fields[10]
                data12 = fields[11]
            result.append(
                Cath_atomic_structure(
                    None,
                    data1,
                    data2,
                    data3,
                    data4,
                    data5,
                    data6,
                    data7,
                    data8,
                    data9,
                    data10,
                    data11,
                    data12))
        file.close()
        os.remove("cath-domain-list.txt")
        return result

    def get_all_structure_x_emd_entry(self):
        emds_id = self.get_all_emds_id('0001', 'inf')
        result = []
        print(len(emds_id))
        k = 0
        for emd in emds_id:
            url = "http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{0}/header/emd-{0}.xml".format(
                emd)
            urllib.request.urlretrieve(url, "{0}.xml".format(emd))
            doc = minidom.parse("{0}.xml".format(emd))

            pdbs = doc.getElementsByTagName("fittedPDBEntryId")
            
            for i in pdbs:
                result.append(Atomic_structure_x_emd_entry(i.firstChild.data, emd))
            os.remove("{0}.xml".format(emd))
            print(str(k) + str(len(result)))
            k+=1
        return result


'''
test = FTP_connection()
test.init_connection()
#directorios = test.get_emds_higher_than_date(datetime.datetime(2019, 6, 1))
directorios = test.get_all_emds_id()
print( len(directorios) )
print(directorios[0])
test.close_connection()
'''
