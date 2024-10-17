query_get_create_partition_table <- function(table_name,partition_by)
{
 query <- paste0('
  CREATE TABLE IF NOT EXISTS public.',table_name,'
  (
      diagnosis_input text COLLATE pg_catalog."default",
      loc_id text COLLATE pg_catalog."default",
      "Country" text COLLATE pg_catalog."default",
      "Type" text COLLATE pg_catalog."default",
      "Year" double precision,
      "Ann. background mortality" double precision,
      "Ann. background population" double precision,
      "Ann. early deaths" double precision,
      "Ann. early deaths (2 diagnosis)" double precision,
      "Ann. early deaths (3 basic care)" double precision,
      "Ann. early deaths (4 best care)" double precision,
      "Ann. early deaths (5 cure)" double precision,
      "Ann. onset deaths" double precision,
      "Ghosts" double precision,
      "Ghosts (delta basic care)" double precision,
      "Ghosts (delta best care)" double precision,
      "Ghosts (delta cure)" double precision,
      "Ghosts (early death)" double precision,
      "Ghosts (onset death)" double precision,
      "Ghosts lever2023" double precision,
      "Ghosts (delta basic care) lever2023" double precision,
      "Ghosts (delta best care) lever2023" double precision,
      "Ghosts (delta cure) lever2023" double precision,
      "Ghosts (early death) lever2023" double precision,
      "Ghosts (onset death) lever2023" double precision,
      "Incidence (1 base)" double precision,
      "Incidence (2 diagnosis)" double precision,
      "Prevalence" double precision,
      "Prevalence stagnant 1970" double precision,
      "Prevalence delay onset by 1 years" double precision,
      "Prevalence delay onset by 3 years" double precision,
      "Prevalence delay onset by 5 years" double precision,
      "Prevalence delay onset by 8 years" double precision,
      "Prevalence delay onset by 13 years" double precision,
      "Ghosts delay onset by 1 years" double precision,
      "Ghosts delay onset by 3 years" double precision,
      "Ghosts delay onset by 5 years" double precision,
      "Ghosts delay onset by 8 years" double precision,
      "Ghosts delay onset by 13 years" double precision,
      age_bracket text COLLATE pg_catalog."default"
    )
    PARTITION BY LIST (',partition_by,');
  ')
 query
}

