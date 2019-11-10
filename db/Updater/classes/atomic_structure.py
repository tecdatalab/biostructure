'''
Created on 22 feb. 2019

@author: luis98
'''
from psycopg2 import sql
import json
import urllib
import os
from generators.LZerD.generate_atom_descriptor import Calculate

class Atomic_structure(object):
    '''
    classdocs
    '''
    __id = None
    __id_code = None
    __parent = None
    __atomic_structure_type_id = None
    __sequence = None
    __atoms = None
    __atoms_count = None
    __aminoacid_count = None
    __png_img_3d = None
    __gif_img_3d = None
    __numbers_descriptor = None
    __domains_strc = []
    __chains_strc = []

    def __init__(
            self,
            id=None,
            id_code=None,
            parent=None,
            atomic_structure_type_id=None,
            sequence=None,
            atoms=None,
            atoms_count=None,
            aminoacid_count=None,
            png_img_3d=None,
            gif_img_3d=None,
            numbers_descriptor=[]):
        self.__id = id
        self.__id_code = id_code
        self.__parent = parent
        self.__atomic_structure_type_id = atomic_structure_type_id
        self.__sequence = sequence
        self.__atoms = atoms
        self.__atoms_count = atoms_count
        self.__aminoacid_count = aminoacid_count
        self.__png_img_3d = png_img_3d
        self.__gif_img_3d = gif_img_3d
        self.__numbers_descriptor = numbers_descriptor
        self.__domains_strc = []
        self.__chains_strc = []
        
    def insert_update_db_complex(self, cur):
        self.insert_update_db(cur)
        for i in self.chains_strc:
            i.parent = self.id
            i.insert_update_db(cur)
        
        for i in self.domains_strc:
            i.parent = self.id
            i.insert_update_db(cur)
        
    def insert_update_db(self, cur):
        cur.execute(
            sql.SQL("SELECT id FROM atomic_structure WHERE id_code = %s;"),[
                self.__id_code])
        result = [record[0] for record in cur]
        if len(result)>0:
            self.__id = result[0]
            self.update_db(cur)
        else:
            self.insert_db(cur)
            
            
    def insert_db(self, cur):
        if self.__id is not None:
            cur.execute(sql.SQL("INSERT INTO atomic_structure(id, id_code, parent, atomic_structure_type_id, sequence, atoms, atoms_count, aminoacid_count, png_img_3d, gif_img_3d, numbers_descriptor) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);"), [
                        self.__id, self.__id_code, self.__parent, self.__atomic_structure_type_id, self.__sequence, self.__atoms, self.__atoms_count, self.__aminoacid_count, self.__png_img_3d, self.__gif_img_3d, json.dumps(self.__numbers_descriptor)])
        else:
            cur.execute(sql.SQL("INSERT INTO atomic_structure(id, id_code, parent, atomic_structure_type_id, sequence, atoms, atoms_count, aminoacid_count, png_img_3d, gif_img_3d, numbers_descriptor) VALUES (default, %s,%s,%s,%s,%s,%s,%s,%s,%s,%s) RETURNING id;"), [
                        self.__id_code, self.__parent, self.__atomic_structure_type_id, self.__sequence, self.__atoms, self.__atoms_count, self.__aminoacid_count, self.__png_img_3d, self.__gif_img_3d, json.dumps(self.__numbers_descriptor)])
            self.__id = [record for record in cur][0]

    
    def update_db(self, cur):
        if self.__id is not None:
            cur.execute(sql.SQL("UPDATE atomic_structure SET parent = %s, atomic_structure_type_id = %s, sequence = %s, atoms = %s, atoms_count = %s, aminoacid_count = %s, png_img_3d = %s, gif_img_3d = %s, numbers_descriptor = %s WHERE id_code = %s;"), [
                        self.__parent, self.__atomic_structure_type_id, self.__sequence, self.__atoms, self.__atoms_count, self.__aminoacid_count, self.__png_img_3d, self.__gif_img_3d, json.dumps(self.__numbers_descriptor),self.__id_code])
        else:
            cur.execute(sql.SQL("UPDATE atomic_structure SET parent = %s, atomic_structure_type_id = %s, sequence = %s, atoms = %s, atoms_count = %s, aminoacid_count = %s, png_img_3d = %s, gif_img_3d = %s, numbers_descriptor = %s WHERE id = %s;"), [
                        self.__parent, self.__atomic_structure_type_id, self.__sequence, self.__atoms, self.__atoms_count, self.__aminoacid_count, self.__png_img_3d, self.__gif_img_3d, json.dumps(self.__numbers_descriptor),self.__id])
    
    def __different_elements(self, list):
        dicc = {}
        count = 0
        for i in list:
            if dicc.get(i) is None:
                count += 1
                dicc[i] = "ok"
        return count
    
    def __aminoacid_list(self, list):
        elements = []
        for i in list:
            temp = i.split()
            temp_add = [temp[3],temp[5]]
            if temp_add not in elements:
                elements.append(temp_add)
        
        result = []
        
        for i in elements:
            result.append(i[0])
        
        return result

    def create_by_online_file(self, id_code, domains):
        
        file_url = "https://files.rcsb.org/download/{0}.pdb".format(str(id_code))
        urllib.request.urlretrieve(file_url, 'atomic-strcuture.txt')
        file = open("atomic-strcuture.txt", "r")

        atom_lines = []
        atom_aminoacid = []
        dic_atom_aminoacid = {}

        for line in file:
            line_split = line.split()
            if line_split[0] == "ATOM":
                atom_lines.append(line)
            elif line_split[0] == "SEQRES":
                atom_aminoacid += line_split[4:-1]
                
                if dic_atom_aminoacid.get(line_split[2]) == None:
                    dic_atom_aminoacid[line_split[2]] = []
                    dic_atom_aminoacid[line_split[2]] += line_split[4:-1]
        
        # Create complex
        self.id_code = id_code
        self.atomic_structure_type_id = 1
        self.sequence = "\n".join(atom_aminoacid)
        self.atoms = "".join(atom_lines)
        self.atoms_count = len(atom_lines)
        self.aminoacid_count = self.__different_elements(atom_aminoacid)
        self.numbers_descriptor = self.calculate_descriptors("".join(atom_lines))

        # Create chains
        actual_chain = ""
        temp_atom_lines = []
        for i in range(len(atom_lines)):
            temp = atom_lines[i].split()[4]
            if (temp != actual_chain and actual_chain != "") or (i+1>=len(atom_lines)):
                temp_add = Atomic_structure(
                    None,
                    str(id_code+actual_chain),
                    id_code,
                    2,
                    "\n".join(dic_atom_aminoacid[actual_chain]),
                    "".join(temp_atom_lines),
                    len(temp_atom_lines),
                    len(dic_atom_aminoacid[actual_chain]),
                    None,
                    None,
                    self.calculate_descriptors("".join(temp_atom_lines)))
                self.__chains_strc.append(temp_add)
                temp_atom_lines = []
            else:
                temp_atom_lines.append(atom_lines[i])
            actual_chain = temp
        # Create domains
        if domains is not None:
            for i in domains:
                add_atoms = []
                for k in i.sequences:
                    add_atoms += atom_lines[k[0] - 1:k[1]]
                    
                atom_aminoacid_domain = self.__aminoacid_list(add_atoms)

                temp = Atomic_structure(
                    None,
                    i.domain,
                    id_code,
                    3,
                    "\n".join(atom_aminoacid_domain),
                    "".join(add_atoms),
                    len(add_atoms),
                    len(atom_aminoacid_domain),
                    None,
                    None,
                    self.calculate_descriptors("".join(add_atoms)))
                self.__domains_strc.append(temp)
        file.close()
        os.remove("atomic-strcuture.txt")

    def calculate_descriptors(self, atoms):
        temp = Calculate()
        return temp.calculate_descriptors(atoms)    
        
    def get_id(self):
        return self.__id

    def get_id_code(self):
        return self.__id_code

    def get_parent(self):
        return self.__parent

    def get_atomic_structure_type_id(self):
        return self.__atomic_structure_type_id

    def get_sequence(self):
        return self.__sequence

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

    def set_id(self, value):
        self.__id = value

    def set_id_code(self, value):
        self.__id_code = value

    def set_parent(self, value):
        self.__parent = value

    def set_atomic_structure_type_id(self, value):
        self.__atomic_structure_type_id = value

    def set_sequence(self, value):
        self.__sequence = value

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

    def del_id(self):
        del self.__id

    def del_id_code(self):
        del self.__id_code

    def del_parent(self):
        del self.__parent

    def del_atomic_structure_type_id(self):
        del self.__atomic_structure_type_id

    def del_sequence(self):
        del self.__sequence

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

    def __str__(self):
        return("id: {0}".format(self.id) + "\n" +
               "id_code: {0}".format(self.id_code) + "\n" +
               "parent: {0}".format(self.parent) + "\n" +
               "atomic_structure_type_id: {0}".format(self.atomic_structure_type_id) + "\n" +
               #"sequence: {0}".format(self.sequence) + "\n" +
               #"atoms: {0}".format(self.atoms) +
               "atoms_count: {0}".format(self.atoms_count) + "\n" +
               "aminoacid_count: {0}".format(self.aminoacid_count) + "\n" +
               "png_img_3d: {0}".format(self.png_img_3d) + "\n" +
               "gif_img_3d: {0}".format(self.gif_img_3d) + "\n" +
               "numbers_descriptor: {0}".format(self.numbers_descriptor))

    id = property(get_id, set_id, del_id, "id's docstring")
    id_code = property(get_id_code, set_id_code, del_id_code, "id_code's docstring")
    parent = property(get_parent, set_parent, del_parent, "parent's docstring")
    atomic_structure_type_id = property(
        get_atomic_structure_type_id,
        set_atomic_structure_type_id,
        del_atomic_structure_type_id,
        "atomic_structure_type_id's docstring")
    sequence = property(
        get_sequence,
        set_sequence,
        del_sequence,
        "sequence's docstring")
    atoms = property(get_atoms, set_atoms, del_atoms, "atoms's docstring")
    atoms_count = property(
        get_atoms_count,
        set_atoms_count,
        del_atoms_count,
        "atoms_count's docstring")
    aminoacid_count = property(
        get_aminoacid_count,
        set_aminoacid_count,
        del_aminoacid_count,
        "aminoacid_count's docstring")
    png_img_3d = property(
        get_png_img_3d,
        set_png_img_3d,
        del_png_img_3d,
        "png_img_3d's docstring")
    gif_img_3d = property(
        get_gif_img_3d,
        set_gif_img_3d,
        del_gif_img_3d,
        "gif_img_3d's docstring")
    numbers_descriptor = property(
        get_numbers_descriptor,
        set_numbers_descriptor,
        del_numbers_descriptor,
        "numbers_descriptor's docstring")
    domains_strc = property(
        get_domains_strc,
        set_domains_strc,
        del_domains_strc,
        "domains_strc's docstring")
    chains_strc = property(
        get_chains_strc,
        set_chains_strc,
        del_chains_strc,
        "chains_strc's docstring")
    numbers_descriptor = property(
        get_numbers_descriptor,
        set_numbers_descriptor,
        del_numbers_descriptor,
        "numbers_descriptor's docstring")
