# 对生存分析结果进行FDR矫正
# usage: Rscript do_adjust.R ${file_path}

suppressMessages({
  library(data.table)
  library(tidyverse)
})

args<-commandArgs(TRUE)

fpath <- args[1]

# fpath <- "/home/wuzj/survival/03.result/methylation/filter_res/ACC_survival_result_filtered"
# fpath <- "/home/wuzj/survival/03.result/protein/L4/ACC_survival_result_with_gene"
# fpath <- "/home/wuzj/survival/03.result/res_DFI/APA/ACC_survival_result"

data <- fread(fpath, header = T, sep = "\t", data.table = F)


if ("multicoxp" %in% colnames(data)) {
  data <- rename(data, "MultiCoxp"="multicoxp")
}


# 不过滤直接矫正
## FDR矫正
data <- data %>%
  mutate(
    KMp.adj = p.adjust(KMp, method = "fdr"),
    UniCoxp.adj = p.adjust(UniCoxp, method = "fdr"),
    MultiCoxp.adj = p.adjust(MultiCoxp, method = "fdr")
  )

## Bonferroni矫正
# data <- data %>%
#   mutate(
#     KMp.adj = p.adjust(KMp, method = "bonferroni"),
#     UniCoxp.adj = p.adjust(UniCoxp, method = "bonferroni"),
#     MultiCoxp.adj = p.adjust(MultiCoxp, method = "bonferroni")
#   )

# if ("pvalHRmedOS" %in% colnames(data)) {
#   data <- data %>%
#     mutate(
#       MultiCoxp.adj = p.adjust(pvalHRmedOS, method = "fdr")
#     )
# } else {
#   warning(paste("Column 'pvalHRmedOS' not found in:", fpath))
# }

output <- paste0(fpath,"_fdr")

fwrite(data, output, sep = "\t", quote = F, row.names = F, col.names = T)


# 过滤再矫正
# if ("aa" %in% colnames(data)) {
#   filtered_data <- data %>%
#     # 任一分组下样本大小低于总数的5%，或高于总数的95%
#     filter(if_else(AA>0,(AA > N*0.05 & AA < N*0.95),T) &
#              if_else(Aa>0,(Aa > N*0.05 & Aa < N*0.95),T) &
#              if_else(aa>0,(aa > N*0.05 & aa < N*0.95),T))
# } else {
#   filtered_data <- data %>%
#     # 任一分组下样本大小低于总数的5%，或高于总数的95%
#     filter((AA > N*0.05 & AA < N*0.95) &
#              (Aa > N*0.05 & Aa < N*0.95))
# }
# 
# 
# filtered_data <- filtered_data %>%
#   mutate(
#     KMp.adj = p.adjust(KMp, method = "fdr"),
#     UniCoxp.adj = p.adjust(UniCoxp, method = "fdr"),
#     MultiCoxp.adj = p.adjust(MultiCoxp, method = "fdr")
#   )
# 
# fwrite(filtered_data, paste0(fpath, "_filtered"), sep = "\t", quote = F, row.names = F, col.names = T)



