# Creates table for Aida's thesis
# Got a reminder on 01/07/2021 to make a table of the biopsies of patients
library("readxl")
library("stringr")
library("purrr")
library("lubridate")
library("dplyr")
library("tidyr")

# RNAseq ####
db <- readxl::read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = "n.a.")
# Using this data as there shouldn't be no problem with w0
conn <- gzfile("data/voom.RNAseq.data.all.cal.noduplications.tsv.gz")
rna <- read.table(conn, sep = "\t", check.names = FALSE)
rna <- rna[ , !grepl(" reseq$", colnames(rna))] # Remove as said
db2 <- db[db$Sample_id %in% colnames(rna), ]
# Clean unnecessary columns
keep_cols <- c(colnames(db2)[1:13], colnames(db2)[20:58])
db2 <- unique(db2[, keep_cols])

#Same metadata than samples ####
stopifnot(nrow(db2) == ncol(rna))
db2 <- db2[match(db2$Sample_id, colnames(rna)), ]

db2 <- mutate(db2, Week = case_when(week == 0 & IBD != "ctrl" ~ "0",
                                    week != 0 & IBD != "ctrl" ~ "14/46",
                                    IBD == "ctrl" ~ "ctrl",
                                    TRUE ~ NA_character_)) %>%
  filter(week %in% c("0", "ctrl"))

db2 %>%
  group_by(IBD, Gender, `Diagnostic age`, sample_location) %>%
  count()


pc <- function(x, column) {
  x %>% group_by(IBD) %>%
    count({{column}}) %>%
    pivot_wider(names_from = {{column}}, values_from = n) %>%
    arrange(forcats::fct_relevel(IBD, "ctrl", "UC", "CD"))
}
# On the study only colon samples were used
db2 %>% filter(sample_location == "colon") %>% pc(Gender)
db2 %>% filter(sample_location == "colon") %>% pc(Ulcers)
db2 %>% filter(sample_location == "colon") %>% pc(`Diagnostic age`)
db2 %>% filter(sample_location == "colon") %>% pc(sample_location)
db2 %>% filter(sample_location == "colon") %>% pc(IBD)
db2 %>%
  filter(sample_location == "colon") %>%
  group_by(IBD) %>%
  summarize(Patients = n_distinct(Pacient_id)) %>%
  arrange(forcats::fct_relevel(IBD, "ctrl", "UC", "CD"))
db2 %>%
  filter(sample_location == "colon") %>%
  group_by(IBD) %>%
  summarize(score = median(CDEIS_partial, na.rm = TRUE))
db2 %>%
  filter(sample_location == "colon") %>%
  group_by(IBD) %>%
  summarize(score = sd(CDEIS_partial, na.rm = TRUE))
db2 %>%
  filter(sample_location == "colon") %>%
  group_by(IBD) %>%
  summarize(score = median(`segment endoscopic MAYO UC`, na.rm = TRUE))
db2 %>%
  filter(sample_location == "colon") %>%
  group_by(IBD) %>%
  summarize(score = sd(`segment endoscopic MAYO UC`, na.rm = TRUE))
