/*
 * Author: Danny Xie Li.
 * Purpose: Creation of the table binnacle
 * to store a binnacle related to updater  
 * emd and pdb.
 */
 
CREATE TABLE segment_entry(
    id_segment SERIAL PRIMARY KEY,
	map_id INT REFERENCES map_information(id),
	algorithm_id INT REFERENCES segment_algorithm(id_algorithm),
	countour_level FLOAT8,
	volume FLOAT8,
	path_map TEXT
);

COMMENT ON TABLE segment_entry IS 'Stores the generate segments by an algorithm of a map.';