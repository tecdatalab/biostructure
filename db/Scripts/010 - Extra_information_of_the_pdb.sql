/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the tables that store the data
 * of the pdb (cath), but that are not regularly consulted.
 *
 */

CREATE TABLE cath_x_atomic_structure(
	id SERIAL PRIMARY KEY,
	cath_domain_name TEXT NOT NULL,
 	class_number INT NOT NULL,
	arquitecture_number INT NOT NULL,
	topology_number INT NOT NULL,
	homologous_superfamiy_number INT NOT NULL,
	S35_sequence_cluster_number INT NOT NULL,
	S60_sequence_cluster_number INT NOT NULL,
	S95_sequence_cluster_number INT NOT NULL,
	S100_sequence_cluster_number INT NOT NULL,
	S100_sequence_count_number INT NOT NULL,
	domain_length INT NOT NULL,
	structure_resolution FLOAT NOT NULL,
    atomic_structure_id INT REFERENCES atomic_structure(id) NOT NULL,
	chain_character TEXT NOT NULL,
	domain_number INT NOT NULL
);

COMMENT ON TABLE atomic_structure IS 'Storage of the cath associated with the structures.';