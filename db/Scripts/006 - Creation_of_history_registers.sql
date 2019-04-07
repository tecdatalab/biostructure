/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of tables for the management of updates.
 *
 */

CREATE TABLE time_stamp(
	emd_entry_id INT REFERENCES emd_entry(id) PRIMARY KEY,
	modification DATE NOT NULL,
	map_file TEXT NOT NULL,
	xml_file TEXT NOT NULL,
	image_file
);

COMMENT ON TABLE time_stamp IS 'Table to know the date of the last EMD update of the database.';

CREATE TABLE update(
	last_update DATE PRIMARY KEY
);

COMMENT ON TABLE update IS 'Table to know the last day in which the updating procedures were executed.';