/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the tables that store
 * the data of the complex, also as well as the 
 * structures related to the pdb.
 */

CREATE TABLE atomic_structure_type(
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL
);

COMMENT ON TABLE atomic_structure_type IS 'Storage of the different types of atomic structures example as the pdb.';

CREATE TABLE atomic_structure(
	id SERIAL PRIMARY KEY,
	id_code TEXT UNIQUE NOT NULL,
	parent INT REFERENCES atomic_structure(id),
    atomic_structure_type_id INT REFERENCES atomic_structure_type(id) NOT NULL,
    sequence TEXT NOT NULL,
    atoms TEXT NOT NULL,
    atoms_count INT NOT NULL,
    aminoacid_count INT NOT NULL,
    png_img_3d TEXT,
    gif_img_3d TEXT,
    numbers_descriptor JSON NOT NULL
);

COMMENT ON TABLE atomic_structure IS 'Storage of the of atomic structures example as the pdb.';

CREATE TABLE atomic_structure_x_emd_entry(
    atomic_structure_id INT REFERENCES atomic_structure(id) NOT NULL,
    emd_entry_id INT REFERENCES emd_entry(id) NOT NULL,
    PRIMARY KEY (atomic_structure_id,emd_entry_id)
);

COMMENT ON TABLE atomic_structure_x_emd_entry IS 'Storage of the connection between emd entry and pdb.';
