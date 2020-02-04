/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the tables that store the data
 * of the atomic structure (cath), but that are not regularly consulted.
 *
 */

CREATE TABLE cath_atomic_structure(
    atomic_structure_id INT REFERENCES atomic_structure(id) PRIMARY KEY,
 	class_number INT NOT NULL,
	architecture_number INT NOT NULL,
	topology_number INT NOT NULL,
	homologous_superfamily_number INT NOT NULL,
	S35_sequence_cluster_number INT NOT NULL,
	S60_sequence_cluster_number INT NOT NULL,
	S95_sequence_cluster_number INT NOT NULL,
	S100_sequence_cluster_number INT NOT NULL,
	S100_sequence_count_number INT NOT NULL,
	domain_length INT NOT NULL,
	structure_resolution FLOAT NOT NULL
);

COMMENT ON TABLE cath_atomic_structure IS 'Storage of the cath associated with the structures.';
