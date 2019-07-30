/*
 * Author: Jose Salas Bonilla.
 * Purpose: Insert default values to representation, user_role and volume_filter tables.
 * Modifier: Luis José Castillo Valverde
 * Purpose: Add initial inserts that were in the updater.
 */

INSERT INTO representation VALUES ('1', 'EMDB contour'), 
('2', 'EMDB contour + 1/3 core'), 
('3', 'EMDB contour + 2/3 core'), 
('4', 'EMDB contour + 1/3 + 2/3 core'), 
('5', 'EMDB contour + 1 std dev');

INSERT INTO user_role VALUES ('1', 'User'), ('2', 'Admin');

INSERT INTO volume_filter VALUES ('0','Off'), ('1','On');

INSERT INTO type_descriptor VALUES (1,'EMDB Contour', 'EMDB Contour'),
(2, 'EMDB Contour + 1/3 core', 'EMDB Contour + 1/3 core'),
(3, 'EMDB Contour + 2/3 core', 'EMDB Contour + 2/3 core'),
(4, 'MDB Contour + 1/3 + 2/3 core', 'EMDB Contour + 1/3 + 2/3 core'),
(5, 'EMDB Contour + 1 std dev', 'EMDB Contour + 1 std dev');
