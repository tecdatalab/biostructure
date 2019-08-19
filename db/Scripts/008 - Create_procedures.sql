/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of procedures for the calculation of the Euclidean distance.
 *
 */

CREATE OR REPLACE FUNCTION euclidean_distance(numbers_emd_1 JSON, numbers_emd_2 JSON, len INTEGER) 
   RETURNS FLOAT AS $$ 
DECLARE
   result FLOAT := 0 ;
   temp_1 FLOAT := 0 ;
   temp_2 FLOAT := 0 ; 
BEGIN

   WHILE len > 0 LOOP
      select (numbers_emd_1->len-1) into temp_1;
      select (numbers_emd_2->len-1) into temp_2;

      result :=result + power((temp_1-temp_2),2); 
      len := len-1;
   END LOOP ; 
   
   RETURN sqrt(result) ;
END;
$$ LANGUAGE plpgsql; 


CREATE OR REPLACE FUNCTION top_distance(emd_id_p INTEGER, type_descriptor INTEGER, top_can INTEGER) 
  RETURNS TABLE (distance   FLOAT
               , emd_id   INT) AS
  $func$
  DECLARE
  len INTEGER := 0;
  numbers_emd JSON;
  BEGIN
    SELECT json_array_length(numbers), numbers INTO len,numbers_emd FROM descriptor WHERE emd_entry_id = emd_id_p AND type_descriptor_id = type_descriptor;
    RETURN QUERY
    SELECT "euclidean_distance_temp"(numbers,numbers_emd,len) as distance, emd_entry_id 
  FROM descriptor
  WHERE type_descriptor_id = type_descriptor AND emd_entry_id != emd_id_p
  ORDER BY distance
  LIMIT top_can;
  END
$func$  LANGUAGE plpgsql;