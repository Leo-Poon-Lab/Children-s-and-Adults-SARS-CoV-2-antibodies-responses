library(tidyverse)
library(readxl)

# data hls
data_raw <- read_xlsx("../data/210916CC ELISA Children Adults MASTER for Haogao_edit.xlsx")
names(data_raw)[1:9] <- c("Group", "ID_group", "ID", "ID_hospital", "Age_group", "Sex", "Age", "Symptoms", "Time")

fill_with_zero <- function(x){
	x[is.na(x)] <- 0
	return(x)
}

data_raw_antibody <- data_raw %>% select(!contains("Avidity"))
data_raw_avidity <- data_raw %>% select(contains("Avidity"))
data_raw_avidity <- data_raw_avidity %>% mutate_all(fill_with_zero)
data_raw_filled <- bind_cols(data_raw_antibody, data_raw_avidity)

write_csv(data_raw_filled, "../data/data_raw_filled.csv")
