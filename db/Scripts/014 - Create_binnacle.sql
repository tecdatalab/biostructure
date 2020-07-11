/*
 * Author: Danny Xie Li.
 * Purpose: Creation of the table binnacle
 * to store a binnacle related to updater  
 * emd and pdb.
 */

CREATE TABLE binnacle_emd(
    binnacle_emd_id SERIAL PRIMARY KEY,
	emd_id INT,
	last_update TIMESTAMP REFERENCES updater(last_update) NOT NULL,
	attempt INT NOT NULL
);

CREATE TABLE binnacle_pdb(
    binnacle_pdb_id SERIAL PRIMARY KEY,
	pdb_id TEXT,
	last_update TIMESTAMP REFERENCES updater(last_update) NOT NULL,
	attempt INT NOT NULL
);

CREATE TABLE updater(
    last_update TIMESTAMP PRIMARY KEY
);

CREATE OR REPLACE PROCEDURE insert_updater(TIMESTAMP)
LANGUAGE plpgsql    
AS $$
BEGIN
    INSERT INTO updater (last_update) 
    VALUES ($1);
END;
$$;

CREATE OR REPLACE PROCEDURE insert_binnacle_emd(INT, TIMESTAMP, INT)
LANGUAGE plpgsql    
AS $$
BEGIN
    INSERT INTO binnacle_emd (emd_id, last_update, attempt) 
    VALUES ($1, $2, $3);
END;
$$;

CREATE OR REPLACE PROCEDURE insert_binnacle_pdb(TEXT, TIMESTAMP, INT)
LANGUAGE plpgsql    
AS $$
BEGIN
    INSERT INTO binnacle_pdb (pdb_id, last_update, attempt) 
    VALUES ($1, $2, $3);
END;
$$;

CREATE OR REPLACE PROCEDURE update_binnacle_emd(INT, TIMESTAMP, INT)
LANGUAGE plpgsql    
AS $$
BEGIN
    UPDATE binnacle_emd   
    SET attempt = (SELECT attempt FROM binnacle_emd WHERE (emd_id = $1 AND last_update = $2)) + $3
    WHERE emd_id = $1 AND last_update = $2;
END;
$$;

CREATE OR REPLACE PROCEDURE update_binnacle_pdb(TEXT, TIMESTAMP, INT)
LANGUAGE plpgsql    
AS $$
BEGIN
    UPDATE binnacle_pdb   
    SET attempt = (SELECT attempt FROM binnacle_pdb WHERE (pdb_id = $1 AND last_update = $2)) + $3
    WHERE pdb_id = $1 AND last_update = $2;
END;
$$;

CREATE OR REPLACE FUNCTION get_attempt_emd(INT, TIMESTAMP) returns INT as $$
DECLARE
    attempt_v INTEGER;
BEGIN
    SELECT attempt INTO attempt_v FROM binnacle_emd WHERE emd_id = $1 AND last_update = $2;
    RETURN attempt_v;
END;
$$ language plpgsql;

CREATE OR REPLACE FUNCTION get_attempt_pdb(TEXT, TIMESTAMP) returns INT as $$
DECLARE
    attempt_v INTEGER;
BEGIN
    SELECT attempt INTO attempt_v FROM binnacle_pdb WHERE pdb_id = $1 AND last_update = $2;
    RETURN attempt_v;
END;
$$ language plpgsql;