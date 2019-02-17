/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the tables associated with the storage of PDBs.
 *
 */

CREATE TABLE pdb_entry(
	id SERIAL PRIMARY KEY,
	pdb TEXT NOT NULL
);

COMMENT ON TABLE pdb_entry IS 'Stores the different PDBs.';

CREATE TABLE pdb_entry_x_emd_entry(
	pdb_entry_id INT REFERENCES pdb_entry(id) NOT NULL,
	emd_entry_id INT REFERENCES emd_entry(id) NOT NULL, 
	PRIMARY KEY (emd_entry_id,pdb_entry_id)
);

COMMENT ON TABLE pdb_entry_x_emd_entry IS 'Associate PDBs with EMDs in a many-to-many relationship.';
