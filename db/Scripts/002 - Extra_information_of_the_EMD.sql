/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the tables that store the data
 * of the EMD, but that are not regularly consulted.
 *
 */

CREATE TABLE map_information(
	id SERIAL PRIMARY KEY,
	file_information JSON NOT NULL,
	data_type TEXT NOT NULL,
	num_columns INT NOT NULL,
	num_rows INT NOT NULL,
	num_sections INT NOT NULL,
	origin_col INT NOT NULL,
	origin_row INT NOT NULL,
	origin_sec INT NOT NULL,
	limit_col INT NOT NULL,
	limit_row INT NOT NULL,
	limit_sec INT NOT NULL,
	spacing_col INT NOT NULL,
	spacing_row INT NOT NULL,
	spacing_sec INT NOT NULL,
	cell_a JSON NOT NULL,
	cell_b JSON NOT NULL,
	cell_c JSON NOT NULL,
	cell_alpha JSON NOT NULL,
	cell_beta JSON NOT NULL,
	cell_gamma JSON NOT NULL,
	axis_order_fast CHAR(1) NOT NULL,
	axis_order_medium CHAR(1) NOT NULL,
	axis_order_slow CHAR(1) NOT NULL,
	minimum FLOAT8 NOT NULL,
	maximum FLOAT8 NOT NULL,
	average FLOAT8 NOT NULL,
	std FLOAT8 NOT NULL,
	space_group_number INT NOT NULL,
	datails TEXT NOT NULL,
	pixel_x FLOAT8 NOT NULL,
	pixel_y FLOAT8 NOT NULL,
	pixel_z FLOAT8 NOT NULL,
	countour_level FLOAT8 NOT NULL,
	annotation_details TEXT
);

COMMENT ON TABLE map_information IS 'Stores the map information of an EMD.';