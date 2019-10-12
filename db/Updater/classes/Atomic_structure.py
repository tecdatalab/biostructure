'''
Created on 22 feb. 2019

@author: luis98
'''
from psycopg2 import sql
import json
import urllib
import os
from email._header_value_parser import Domain


class Atomic_structure(object):
    '''
    classdocs
    '''
    __id = None
    __pdb = None
    __parent = None
    __atomic_structure_type_id = None
    __secuence = None
    __atoms = None
    __atoms_count = None
    __aminoacid_count = None
    __png_img_3d = None
    __gif_img_3d = None
    __numbers_descriptor = None
    __domains_strc = []
    __chains_strc = []

    def __init__(self, id = None, pdb = None, parent = None, atomic_structure_type_id = None, secuence = None, atoms = None, atoms_count = None, aminoacid_count = None, png_img_3d = None, gif_img_3d = None, numbers_descriptor = [], domains_strc = [], chains_strc = []):
        self.__id = id
        self.__pdb = pdb
        self.__parent = parent
        self.__atomic_structure_type_id = atomic_structure_type_id
        self.__secuence = secuence
        self.__atoms = atoms
        self.__atoms_count = atoms_count
        self.__aminoacid_count = aminoacid_count
        self.__png_img_3d = png_img_3d
        self.__gif_img_3d = gif_img_3d
        self.__numbers_descriptor = numbers_descriptor
        self.__domains_strc = domains_strc
        self.__chains_strc = chains_strc
        
    def insert_db(self, cur):
        if self.__id != None:
            cur.execute(sql.SQL("INSERT INTO atomic_structure(id, pdb, parent, atomic_structure_type_id, secuence, atoms, atoms_count, aminoacid_count, png_img_3d, gif_img_3d, numbers_descriptor) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);")
            , [self.__id, self.__pdb, self.__parent, self.__atomic_structure_type_id, self.__secuence, self.__atoms, self.__atoms_count, self.__aminoacid_count, self.__png_img_3d, self.__gif_img_3d, json.dumps(self.__numbers)])
        else:
            cur.execute(sql.SQL("INSERT INTO atomic_structure(default, pdb, parent, atomic_structure_type_id, secuence, atoms, atoms_count, aminoacid_count, png_img_3d, gif_img_3d, numbers_descriptor) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s) RETURNING id;")
            , [self.__pdb, self.__parent, self.__atomic_structure_type_id, self.__secuence, self.__atoms, self.__atoms_count, self.__aminoacid_count, self.__png_img_3d, self.__gif_img_3d, json.dumps(self.__numbers)])
            self.__id = [record for record in cur][0]
    
    def __different_elements(self, list):
        dicc = {}
        count = 0
        for i in list:
            if dicc.get(i) == None:
                count += 1
                dicc[i] = "ok"
        return count
    
    def create_by_online_file(self, pdb, domains):
        file_url = "https://files.rcsb.org/download/{0}.pdb".format(str(pdb))
        urllib.request.urlretrieve(file_url, 'atomic-strcuture.txt')
        file = open("atomic-strcuture.txt", "r")
        
        atom_lines = []
        atom_aminoacid = []
        
        for line in file:
            line_split = line.split()
            if line_split[0] == "ATOM":
                atom_lines.append(line)
            elif line_split[0] == "SEQRES":
                atom_aminoacid+=line_split[4:-1]
        
        # Create complex
        self.pdb = pdb
        self.atomic_structure_type_id = 1
        self.secuence = "\n".join(atom_aminoacid)
        self.atoms = "".join(atom_lines)
        self.atoms_count = len(atom_lines)
        self.aminoacid_count = self.__different_elements(atom_aminoacid)
        self.numbers_descriptor = []

        # Create chains
        actual_chain = ""
        temp_atom_lines = []
        for i in range(len(atom_lines)):
            temp = atom_lines[i].split()[4]
            if temp!= actual_chain and actual_chain!="":
                temp_add = Atomic_structure(
                    None,
                    temp,
                    None,
                    2,
                    "\n".join([]),
                    "".join(temp_atom_lines),
                    len(temp_atom_lines),
                    self.__different_elements(temp_atom_lines),
                    None,
                    None,
                    [],
                    [],
                    [])
                self.chains_strc.append(temp)
                temp_atom_lines = []
            else:
                temp_atom_lines.append(atom_lines[i])
            actual_chain = temp
        # Create domains
        if domains != None :
                for i in domains:
                    add_secuence = []
                    add_atoms = []
                    for k in i.sequences:
                        add_secuence += atom_aminoacid[k[0] - 1:k[1]]
                        add_atoms += atom_lines[k[0] - 1:k[1]]
                        
                    temp = Atomic_structure(
                        None,
                        i.domain,
                        None,
                        3,
                        "\n".join(add_secuence),
                        "".join(add_atoms),
                        len(add_atoms),
                        self.__different_elements(add_secuence),
                        None,
                        None,
                        [],
                        [],
                        [])
                    self.domains_strc.append(temp)
        file.close()
        os.remove("atomic-strcuture.txt")
        
    def get_id(self):
        return self.__id

    def get_pdb(self):
        return self.__pdb

    def get_parent(self):
        return self.__parent

    def get_atomic_structure_type_id(self):
        return self.__atomic_structure_type_id

    def get_secuence(self):
        return self.__secuence

    def get_atoms(self):
        return self.__atoms

    def get_atoms_count(self):
        return self.__atoms_count

    def get_aminoacid_count(self):
        return self.__aminoacid_count

    def get_png_img_3d(self):
        return self.__png_img_3d

    def get_gif_img_3d(self):
        return self.__gif_img_3d

    def get_numbers_descriptor(self):
        return self.__numbers_descriptor

    def get_domains_strc(self):
        return self.__domains_strc

    def get_chains_strc(self):
        return self.__chains_strc

    def get_numbers(self):
        return self.__numbers

    def set_id(self, value):
        self.__id = value

    def set_pdb(self, value):
        self.__pdb = value

    def set_parent(self, value):
        self.__parent = value

    def set_atomic_structure_type_id(self, value):
        self.__atomic_structure_type_id = value

    def set_secuence(self, value):
        self.__secuence = value

    def set_atoms(self, value):
        self.__atoms = value

    def set_atoms_count(self, value):
        self.__atoms_count = value

    def set_aminoacid_count(self, value):
        self.__aminoacid_count = value

    def set_png_img_3d(self, value):
        self.__png_img_3d = value

    def set_gif_img_3d(self, value):
        self.__gif_img_3d = value

    def set_numbers_descriptor(self, value):
        self.__numbers_descriptor = value

    def set_domains_strc(self, value):
        self.__domains_strc = value

    def set_chains_strc(self, value):
        self.__chains_strc = value

    def set_numbers(self, value):
        self.__numbers = value

    def del_id(self):
        del self.__id

    def del_pdb(self):
        del self.__pdb

    def del_parent(self):
        del self.__parent

    def del_atomic_structure_type_id(self):
        del self.__atomic_structure_type_id

    def del_secuence(self):
        del self.__secuence

    def del_atoms(self):
        del self.__atoms

    def del_atoms_count(self):
        del self.__atoms_count

    def del_aminoacid_count(self):
        del self.__aminoacid_count

    def del_png_img_3d(self):
        del self.__png_img_3d

    def del_gif_img_3d(self):
        del self.__gif_img_3d

    def del_numbers_descriptor(self):
        del self.__numbers_descriptor

    def del_domains_strc(self):
        del self.__domains_strc

    def del_chains_strc(self):
        del self.__chains_strc

    def del_numbers(self):
        del self.__numbers
        
    def __str__( self ):
        return("id: {0}".format(self.id) +"\n"+
               "pdb: {0}".format(self.pdb) +"\n"+
               "parent: {0}".format(self.parent) +"\n"+
               "atomic_structure_type_id: {0}".format(self.atomic_structure_type_id) +"\n"+
               "secuence: {0}".format(self.secuence) +"\n"+
               "atoms: {0}".format(self.atoms) +
               "atoms_count: {0}".format(self.atoms_count) +"\n"+
               "aminoacid_count: {0}".format(self.aminoacid_count) +"\n"+
               "png_img_3d: {0}".format(self.png_img_3d) +"\n"+
               "gif_img_3d: {0}".format(self.gif_img_3d) +"\n"+
               "numbers_descriptor: {0}".format(self.numbers_descriptor))

    id = property(get_id, set_id, del_id, "id's docstring")
    pdb = property(get_pdb, set_pdb, del_pdb, "pdb's docstring")
    parent = property(get_parent, set_parent, del_parent, "parent's docstring")
    atomic_structure_type_id = property(get_atomic_structure_type_id, set_atomic_structure_type_id, del_atomic_structure_type_id, "atomic_structure_type_id's docstring")
    secuence = property(get_secuence, set_secuence, del_secuence, "secuence's docstring")
    atoms = property(get_atoms, set_atoms, del_atoms, "atoms's docstring")
    atoms_count = property(get_atoms_count, set_atoms_count, del_atoms_count, "atoms_count's docstring")
    aminoacid_count = property(get_aminoacid_count, set_aminoacid_count, del_aminoacid_count, "aminoacid_count's docstring")
    png_img_3d = property(get_png_img_3d, set_png_img_3d, del_png_img_3d, "png_img_3d's docstring")
    gif_img_3d = property(get_gif_img_3d, set_gif_img_3d, del_gif_img_3d, "gif_img_3d's docstring")
    numbers_descriptor = property(get_numbers_descriptor, set_numbers_descriptor, del_numbers_descriptor, "numbers_descriptor's docstring")
    domains_strc = property(get_domains_strc, set_domains_strc, del_domains_strc, "domains_strc's docstring")
    chains_strc = property(get_chains_strc, set_chains_strc, del_chains_strc, "chains_strc's docstring")
    numbers = property(get_numbers, set_numbers, del_numbers, "numbers's docstring")
 
