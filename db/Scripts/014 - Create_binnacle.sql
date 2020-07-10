CREATE TABLE binnacle_emd(
    emd_id INT PRIMARY KEY,
	last_update DATE NOT NULL,
	attempt INT NOT NULL
);

CREATE TABLE binnacle_pdb(
    pdb_id TEXT PRIMARY KEY,
	last_update DATE NOT NULL,
	attempt INT NOT NULL
);

CREATE TABLE updater(
    last_update DATE PRIMARY KEY
);

CREATE OR REPLACE PROCEDURE update_binnacle_pdb(TEXT, DATE, INT)
LANGUAGE plpgsql    
AS $$
BEGIN
    INSERT INTO binnacle_pdb (pdb_id, last_update, attempt) 
    VALUES ($1, $2, $3)
    ON CONFLICT (pdb_id) DO UPDATE 
    SET last_update = excluded.last_update, 
        attempt = 
        (
        SELECT attempt FROM binnacle_pdb WHERE pdb_id = $1
    ) + 1;
    /*COMMIT*/
END;
$$;

CREATE OR REPLACE PROCEDURE update_binnacle_emd(INT, DATE, INT)
LANGUAGE plpgsql    
AS $$
BEGIN
    INSERT INTO binnacle_emd (emd_id,last_update, attempt) 
    VALUES ($1, $2, $3)
    ON CONFLICT (emd_id) DO UPDATE 
    SET last_update = excluded.last_update, 
        attempt = (SELECT attempt FROM binnacle_emd WHERE emd_id = $1) + 1;
    /*COMMIT*/
END;
$$;

CREATE OR REPLACE FUNCTION get_attempt_pdb(TEXT) returns INT as $$
DECLARE
    attempt_v INTEGER;
BEGIN
    SELECT attempt INTO attempt_v FROM binnacle_pdb WHERE pdb_id = $1;
    RETURN attempt_v;
END;
$$ language plpgsql;

CREATE OR REPLACE FUNCTION get_attempt_emd(INT) returns INT as $$
DECLARE
    attempt_v INTEGER;
BEGIN
    SELECT attempt INTO attempt_v FROM binnacle_emd WHERE emd_id = $1;
    RETURN attempt_v;
END;
$$ language plpgsql;