/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of the database.
 *
 */

CREATE DATABASE biomolecules_db
  WITH OWNER = postgres
       ENCODING = 'UTF8'
       TABLESPACE = pg_default
       LC_COLLATE = 'es_CR.UTF-8'
       LC_CTYPE = 'es_CR.UTF-8'
       CONNECTION LIMIT = -1;