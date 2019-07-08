/*
 * Author: Luis Castillo-Valverde.
 * Purpose: Creation of procedures for the calculation of the Euclidean distance.
 *
 */

 CREATE OR REPLACE FUNCTION euclidean_distance(emd_id_1 INTEGER,emd_id_2 INTEGER, type_descriptor INTEGER) 
   RETURNS FLOAT AS $$ 
DECLARE
   len INTEGER := 0 ; 
   result FLOAT := 0 ;
   temp_1 FLOAT := 0 ;
   temp_2 FLOAT := 0 ; 
BEGIN

   select json_array_length(numbers) into len from descriptor WHERE emd_entry_id = emd_id_1 and type_descriptor_id = type_descriptor;
   
   WHILE len > 0 LOOP
      select (numbers->len-1) into temp_1 from descriptor WHERE emd_entry_id = emd_id_1 and type_descriptor_id = type_descriptor;
      select (numbers->len-1) into temp_2 from descriptor WHERE emd_entry_id = emd_id_2 and type_descriptor_id = type_descriptor;

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
  BEGIN
    RETURN QUERY
    SELECT "euclidean_distance"(emd_entry_id,emd_id_p,type_descriptor) as distance,emd_entry_id 
    from descriptor 
    where type_descriptor_id = type_descriptor and emd_entry_id !=emd_id_p
    order by distance
    limit top_can;
  END
$func$  LANGUAGE plpgsql;
