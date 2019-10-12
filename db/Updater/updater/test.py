'''
Created on 5 oct. 2019

@author: luis98
'''

from connections.FTP_connection import FTP_connection
from classes.Atomic_structure import Atomic_structure


def insert_all_pdb():
    conec_ftp = FTP_connection()
    all_pdb = conec_ftp.get_all_pdb()
    domain_dic = conec_ftp.get_all_cath_domain_boundarie_dic()
    for i in all_pdb:
        temp = Atomic_structure()
        temp.create_by_online_file(i.data1, domain_dic.get(i.data1))
        print(temp.__str__())
        print("\n\n------------Chains------------\n\n")
        for i in temp.chains_strc:
            print(i.__str__())
            print("----------------------------------------------")
        print("\n\n------------Domains------------\n\n")
        for i in temp.domains_strc:
            print(i.__str__())
            print("----------------------------------------------")
        break


# insert_all_pdb()
conec_ftp = FTP_connection()
conec_ftp.init_connection()
result = conec_ftp.get_all_pdb_entry_x_emd_entry(1, 3)
#result = conec_ftp.get_all_cath_chain()
#result = conec_ftp.get_all_cath_domain_boundarie()
#result = conec_ftp.get_all_cath_domain()
#result = conec_ftp.get_all_cath_domain_boundarie_dic()
