#' @export
dbGetSpectrum <- function(db, Exp.ID) {
  stopifnot(length(Exp.ID) == 1)
  
  # Selection of Spectrum by Exp.ID
  dbSpectrum <- db %>% # Selection of database
    tbl(., "fragment") %>% # Selection of Fragment-Table (id, mz, int, experiment_id)
    filter(., experiment_id == Exp.ID) %>% # Selection of one Experiment
    select(., mz, int) %>% # mz+int- Table
    dplyr::collect()
  
  mz_i <- db %>%
    tbl("experiment") %>%
    filter(experiment_id == Exp.ID) %>%
    dplyr::collect() %>%
    .$mz
  
  meta_data <- id_to_name(Exp.ID, db)
  
  attr(dbSpectrum, "precursor_mz") <- mz_i
  attr(dbSpectrum, "comp_name") <- meta_data$name
  attr(dbSpectrum, "cas") <- meta_data$cas
  attr(dbSpectrum, "db_exp_ID") <- Exp.ID
  return(dbSpectrum)
}
id_to_name <- function(id, db) {
  compTbl <- tbl(db, "compound")
  name <- tbl(db, "experiment") %>%
    filter(experiment_id == id) %>%
    inner_join(compTbl, by = "compound_id") %>%
    select(name, cas) %>%
    dplyr::collect()
  name$exp_ID <- id
  name
}