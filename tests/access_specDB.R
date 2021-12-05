# Script to manually connect to Spec DB.


library(ntsworkflow)
library(dplyr)

dbi <- DBI::dbConnect(RSQLite::SQLite(), "~/sqlite_local/MS2_db_v7.db")

compTable <- tbl(dbi, "compound")

compLink <- tbl(dbi, "compGroupComp")

compGroup <- tbl(dbi, "compoundGroup")

all_comps <- compTable %>% left_join(compLink, by = "compound_id") %>% 
  left_join(compGroup, by = "compoundGroup_id") %>% collect()

psm <- subset(all_comps, name.y %in% c("Herbicide", "Insecticide", "Fungicide", "Pesticide"))


pcp <- subset(all_comps, name.y == "Personal_care_product")

ip <-  subset(all_comps, name.y == "Industrial_process")

pharma <-  subset(all_comps, name.y == "Pharmaceutical")

food <- subset(all_comps, name.y == "Food_additive")
