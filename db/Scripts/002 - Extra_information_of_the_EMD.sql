/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the tables that store the data
 * of the EMD, but that are not regularly consulted.
 *
 */

CREATE TABLE map_information(
	id SERIAL PRIMARY KEY,
	file_information JSON,
	data_type TEXT,
	num_columns INT,
	num_rows INT,
	num_sections INT,
	origin_col INT,
	origin_row INT,
	origin_sec INT,
	limit_col INT,
	limit_row INT,
	limit_sec INT,
	spacing_col INT,
	spacing_row INT,
	spacing_sec INT,
	cell_a JSON,
	cell_b JSON,
	cell_c JSON,
	cell_alpha JSON,
	cell_beta JSON,
	cell_gamma JSON,
	axis_order_fast CHAR(1),
	axis_order_medium CHAR(1),
	axis_order_slow CHAR(1),
	minimum FLOAT8,
	maximum FLOAT8,
	average FLOAT8,
	std FLOAT8,
	space_group_number INT,
	details TEXT,
	pixel_x JSON,
	pixel_y JSON,
	pixel_z JSON,
	countour_level FLOAT8,
	annotation_details TEXT,
	volume TEXT[]
);

COMMENT ON TABLE map_information IS 'Stores the map information of an EMD.';
