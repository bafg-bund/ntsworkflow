blah <- DBI::dbConnect(RSQLite::SQLite(), "~/sqlite_local/MS2_db_v9.db")
DBI::dbListTables(blah)
DBI::dbReadTable(blah, "experimentGroup")$name
DBI::dbReadTable(blah, "experimentGroup")
head(DBI::dbReadTable(blah, "retention_time"))
head(DBI::dbReadTable(blah, "expGroupexp"))
head(DBI::dbReadTable(blah, "parameter"))

unique(DBI::dbReadTable(blah, "retention_time")$chrom_method)

rtt <- DBI::dbReadTable(blah, "retention_time")

badComps <- rtt[is.na(rtt$chrom_method), "compound_id"]

tbl(blah, "compound") %>% filter(compound_id %in% badComps) %>% left_join(tbl(blah, "compGroupComp"), by = "compound_id") %>% 
  left_join(tbl(blah, "compoundGroup"), by = "compoundGroup_id") %>% collect() %>% View()

iwant <- c("BfG", "LfU")
iwant <- "LfU"
library(dplyr)
lfuTest <- tbl(blah, "experimentGroup") %>% filter(name %in% iwant) %>% 
  left_join(tbl(blah, "expGroupExp"), by = "experimentGroup_id") %>% 
  select(experiment_id) %>% left_join(tbl(blah, "experiment"), by = "experiment_id") %>% 
  left_join(tbl(blah, "parameter"), by = "parameter_id") %>% 
  left_join(tbl(blah, "compound"), by = "compound_id") %>%
  left_join(tbl(blah, "retention_time"), by = "compound_id") %>% 
  filter(name == "Carbamazepine") %>% 
  collect()

DBI::dbDisconnect(blah)


# need to add predicted retention times to the spectral database
# first view relationship between retentions times from the three different methods
db <- DBI::dbConnect(RSQLite::SQLite(), "~/sqlite_local/MS2_db_v9.db")

rtt <- DBI::dbReadTable(db, "retention_time")
head(rtt)
table(rtt$chrom_method)

rt2 <- tbl(db, "compound") %>% left_join(tbl(db, "retention_time"), by = "compound_id") %>% 
  filter(chrom_method %in% c("dx.doi.org/10.1016/j.chroma.2015.11.014", "lfu_nts_rp1")) %>% 
  select(name, chrom_method, rt) %>% collect()
head(rt2)
rt2[rt2$chrom_method == "dx.doi.org/10.1016/j.chroma.2015.11.014", "chrom_method"] <- "bfg_nts_rp1" 
rt3 <- tidyr::pivot_wider(rt2, names_from = chrom_method, values_from = rt)
rt4 <- subset(rt3, !is.na(lfu_nts_rp1))
rt4 <- subset(rt4, !is.na(bfg_nts_rp1))
library(ggplot2)
library(ggrepel)
mod <- lm(lfu_nts_rp1~bfg_nts_rp1, rt4)
mod
length(resid(mod))
data.frame(row_numresid(mod))
rt4$resid <- resid(mod)

rt5 <- subset(rt3, is.na(bfg_nts_rp1))

ggplot(rt4, aes(bfg_nts_rp1, lfu_nts_rp1, group = name, label = name)) + geom_point(alpha = 0.3, shape = 1) + 
  geom_text_repel(alpha = 0.5, size = 2, max.overlaps = 20) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.5, color = "red") +
  geom_abline(slope = 0.9963, intercept = -0.2, alpha = 0.5, color = "blue") +
  geom_hline(data = rt5, mapping = aes(yintercept = lfu_nts_rp1), alpha = 0.2, size = .1)

# add rows where bfg retention time is not known

# which compounds have only lfu retention time?
rt6 <- tbl(db, "compound") %>% left_join(tbl(db, "retention_time"), by = "compound_id") %>% 
  filter(chrom_method %in% c("dx.doi.org/10.1016/j.chroma.2015.11.014", "lfu_nts_rp1")) %>% 
  select(compound_id, name, chrom_method, rt) %>% collect()
rt6 <- tidyr::pivot_wider(rt6, names_from = chrom_method, values_from = rt)
rt6 <- subset(rt6, is.na(`dx.doi.org/10.1016/j.chroma.2015.11.014`))
max(rtt$ret_time_id)

data.frame(ret_time_id = seq_len(nrow(rt6)) + 1862, rt = rt6$lfu_nts_rp1,
           chrom_method = "dx.doi.org/10.1016/j.chroma.2015.11.014", compound_id = rt6$compound_id, predicted = "TRUE")

noId <- data.frame(rt = rt6$lfu_nts_rp1,
           chrom_method = "dx.doi.org/10.1016/j.chroma.2015.11.014", compound_id = rt6$compound_id, predicted = "TRUE")

DBI::dbAppendTable(db, "retention_time", noId)

DBI::dbDisconnect(db)
