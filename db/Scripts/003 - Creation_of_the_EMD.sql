/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the tables that store
 * the data of the EMD, also the different types of
 * descriptors that an EMD can have are stored.
 *
 */

CREATE TABLE type_descriptor(
	id SERIAL PRIMARY KEY,
	name TEXT UNIQUE NOT NULL,
	description TEXT NOT NULL
);

COMMENT ON TABLE type_descriptor IS 'Stores the different types of existing descriptors.';

CREATE TABLE emd_entry(
	id INT PRIMARY KEY, /* The id corresponds to the 4 digits of EMDB */
	full_name TEXT NOT NULL,
	acronym TEXT NOT NULL,
	volume FLOAT8 NOT NULL,
	resolution FLOAT8 NOT NULL,
	image_url TEXT NOT NULL,
	xml_url TEXT NOT NULL,
	map_url TEXT NOT NULL,
	map_information_id INT REFERENCES map_information(id) NOT NULL
);

COMMENT ON TABLE emd_entry IS 'Stores existing EMDs.';

CREATE TABLE descriptor(
	emd_entry_id INT REFERENCES emd_entry(id) NOT NULL, 
	type_descriptor_id INT REFERENCES type_descriptor(id) NOT NULL, 
	numbers JSON NOT NULL,
	PRIMARY KEY (emd_entry_id,type_descriptor_id)
);

COMMENT ON TABLE descriptor IS 'It relates the EMD with the descriptors, in a one-to-one relationship.';