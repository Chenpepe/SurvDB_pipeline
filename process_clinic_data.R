# 临床数据处理

library(data.table)
library(tidyverse)

fpath.clinic <- "/home/wuzj/survival/01.data/clinic/TCGA-CDR.tsv"

data <- fread(fpath.clinic, header = T, sep = "\t") %>%
  select(c(1:6,25:32)) %>% rename(
  sampleID=bcr_patient_barcode,
  age=age_at_initial_pathologic_diagnosis,
  stage=ajcc_pathologic_tumor_stage
)

## 处理肿瘤分期：按大类转为对应数值
data <- data %>%
  mutate(stage = case_when(
    stage == "Stage 0" ~ 0,
    stage %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC", "IS") ~ 1,
    stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC", "Stage IIA", "Stage IIB", "Stage IIC") ~ 2,
    stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ 3,
    stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ 4,
    TRUE ~ NA
  ))

## 按癌症拆分
split_data <- split(data, data$type)

## 保存
for (type in names(split_data)) {
  fname <- paste0("/home/wuzj/survival/01.data/clinic/", type, ".clinic.tsv")
  fwrite(split_data[[type]], file = fname, sep = "\t", row.names = F)
}






