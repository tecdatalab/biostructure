/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the tables that store the data for each
 * of the queries that are made on the 
 * page: http://kiharalab.org/em-surfer/submit.php.
 *
 */

CREATE TABLE representation(
	id SERIAL PRIMARY KEY,
	name TEXT UNIQUE NOT NULL
);

COMMENT ON TABLE representation IS 'Stores the types of representations of the contour shape.';

CREATE TABLE volume_filter(
	id SERIAL PRIMARY KEY,
	name TEXT UNIQUE NOT NULL
);

COMMENT ON TABLE volume_filter IS 'Stores the states that the volume filter can be (off, on).';

CREATE TABLE search_history(
	id SERIAL PRIMARY KEY,
	date_time TIMESTAMP NOT NULL,
	ip TEXT NOT NULL,
	emd_entry_id INT REFERENCES emd_entry(id) NOT NULL,
	name_file TEXT,
	counter_level FLOAT8,
	representation_id INT REFERENCES representation(id) NOT NULL,
	volume_filter_id INT REFERENCES volume_filter(id) NOT NULL,
	resolution_filter_min FLOAT8,
	resolution_filter_max FLOAT8
);

COMMENT ON TABLE search_history IS 'Stores the significant data of the queries.';

CREATE TABLE benchmark_history(
    id SERIAL PRIMARY KEY,
    date_time TIMESTAMP NOT NULL,
    ip TEXT NOT NULL,
    user_id INT,
    representation_id INT REFERENCES representation(id) NOT NULL,
    volume_filter_id INT REFERENCES volume_filter(id) NOT NULL,
    top_results INT NOT NULL,
    emd_list JSON NOT NULL
);

COMMENT ON TABLE benchmark_history IS 'Stores the significant data of the benchmark queries.';
