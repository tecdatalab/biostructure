/*
 * Author: Danny Xie Li.
 * Purpose: Creation of the table segment_algorithm
 * to store a the algorithms used in Updater for 
 * segmentation.
 */

CREATE TABLE segment_algorithm(
    id_algorithm SERIAL PRIMARY KEY,
	algorithm_name TEXT
);

COMMENT ON TABLE segment_algorithm IS 'Stores an enum of the algorithm used to segment the map.';

/* Insert default values */
INSERT INTO segment_algorithm(algorithm_name) VALUES('default'); COMMIT;