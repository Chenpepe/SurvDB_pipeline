# 数据过滤

#### 导包 ####
suppressMessages({
  library(data.table)
  library(tidyverse)
  # library(R.utils)
})

#### 传入参数 ####
args<-commandArgs(TRUE)
### 癌症类型
cancer_type <- args[1]
### 标志物类型
marker_type <- args[2]
### 结局类型
survival_type <- args[3]
### 标志物数据路径
fpath.marker <- args[4]
### 临床数据路径
fpath.clinic <- args[5]
### 输出路径
output_dir <- args[6]

#### 定义函数 ####
# 计算缺失率
calculate_missing_rate <- function(line) {
  total <- length(line)
  num_missing <- sum(is.na(line))
  missing_rate <- num_missing / total
  return(missing_rate)
}

# calculate_missing_rate(line)

# 计算 MAF


# calculate_maf(line)

# 表型数据列名处理
process_marker_data <- function(data.marker, marker_type) {
  data.marker <- switch(marker_type,
                        "AS" = {
                          data.marker %>%
                            mutate(marker = paste(symbol, splice_type, as_id, sep = "_")) %>%
                            select(-c(1:10)) %>%
                            mutate_all(~na_if(as.character(.), "null"))
                        },
                        "APA" = {
                          data.marker %>%
                            rename(marker = event_id) %>%
                            select(-c(2,3))
                        },
                        "gene" = {
                          data.marker %>%
                            rename(marker = gene)
                        },
                        "mRNA" = {
                          data.marker %>%
                            rename(marker = gene)
                        },
                        "lncRNA" = {
                          data.marker %>%
                            rename(marker = gene)
                        },
                        "miRNA" = {
                          data.marker %>%
                            rename(marker = sample)
                        },
                        "methylation" = {
                          data.marker %>%
                            rename(marker = "sample")
                        },
                        "SNP" = {
                          data.marker %>%
                            # mutate(ID = sub(":.*", "", ID)) %>%
                            rename(marker = ID) %>%
                            rename_with(~sub("-01$", "", .), starts_with("TCGA"))
                        },
                        "CNV" = {
                          data.marker %>%
                            rename(marker = "Gene Symbol")
                        },
                        "protein" = {
                          data.marker %>%
                            rename(marker = "Sample_ID")
                        }
  )
  return(data.marker)
}

# 合并表型数据和临床数据
process_survival_data <- function(marker_type,survival_type, line, data.clinic) {
  ### 生成生存分析数据
  ## 测试
  # line <- data.marker[2,]
  
  ## 宽格式转长格式，用于与临床数据合并
  tmp <- gather(line, key = sampleID)[-1,]
  
  ## sampleID 处理
  tmp <- tmp %>%
    # mutate(sampleID = ifelse(grepl("^TCGA-[A-Z0-9]+-[A-Z0-9]+$", sampleID),
    #                          sampleID, sub("-[^-]*$", "", sampleID))) %>%
    mutate(sampleID = sub("^([^\\-]+-[^\\-]+-[^\\-]+).*", "\\1", sampleID)) %>%
    mutate(sampleID = gsub("_", "-", sampleID),
           value = as.numeric(value))
  
  ## 合并表型数据与临床数据
  merge_data <- merge(tmp, data.clinic, by="sampleID") %>%
    # 去除标志物表型数据缺失的样本记录
    filter(!is.na(value))
  
  ## 分组
  if (marker_type == "SNP") {
    merge_data <- mutate(merge_data, group = value)
  } else if (marker_type == "CNV") {
    merge_data <- mutate(merge_data, group = ifelse(value == 0, "normal",
                                                    ifelse(value > 0, "gain", "loss")))
  } else {
    merge_data <- mutate(merge_data, group = ifelse(value > median(merge_data$value), "high", "low"))
  }
  
  ## 肿瘤分期转换数值
  # merge_data <- mutate(merge_data, stage = if_else(stage == "not reported", NA, stage))
  # 
  # stage_table <- data.frame(
  #   stage = c("stage i", "stage ia", "stage ib", "stage ii", "stage iia", "stage iib", "stage iii", "stage iiia", "stage iiib", "stage iv"),
  #   stage_val = 1:10
  # )
  # merge_data <- merge(merge_data, stage_table, by = "stage", all.x = TRUE)
  # merge_data$stage <- merge_data$stage_val
  
  ## 动态选择生存结局类型
  time_col <- paste0(survival_type, ".time")  # 动态生成生存时间列名
  status_col <- survival_type  # 动态生成生存状态列名
  
  survival_data <- merge_data %>%
    # mutate(status = ifelse(.data[[status_col]] == "alive", 0, 1)) %>% # 动态生成状态列
    mutate(status = .data[[status_col]]) %>%
    ## 去除基因型和临床数据缺失的记录
    filter(!(is.na(value) | is.na(.data[[time_col]]) | is.na(.data[[status_col]]))) %>%
    ## 将生存时间在10年截断
    mutate(status = ifelse(.data[[time_col]] > 3650, 0, status)) %>%
    mutate(time = ifelse(.data[[time_col]] > 3650, 3650, .data[[time_col]])) %>%
    ## 将时间单位改为年
    mutate(time = time / 365)
  
  return(survival_data)
}

# 数据过滤
filter_data <- function(marker_type, line, data.clinic) {
  # line <- data.marker[1,]
  marker <- line$marker
  reasons <- c()
  
  ## 根据数据量大小，使用不同缺失阈值
  missing_rate <- calculate_missing_rate(line)
  
  if (marker_type %in% c("APA","AS","gene","mRNA","lncRNA","miRNA","CNV","protein")) {
    cutoff_missing_date <- 0.1
  } else if (marker_type %in% c("SNP","methylation")) {
    cutoff_missing_date <- 0.05
  }
  
  if (missing_rate > cutoff_missing_date) {
    reasons <- c(reasons, "missing_rate")
  }
  
  ## 合并临床数据
  survival_data <- process_survival_data(marker_type,survival_type, line, data.clinic)
  
  
  ## 样本数量统计（step2）
  num_sample <- dim(survival_data)[1]
  if (num_sample != 0) {
    num_group <- length(unique(survival_data$group))
    min_num_group <- min(table(survival_data$group))
    num_surv <- sum(survival_data$status == 1, na.rm = TRUE)

    ## 统计意义
    if (num_sample < 30) { # 有效样本量 < 30
      reasons <- c(reasons, "num_sample")
    }
    if (num_group < 2) { # 分组数 < 2
      reasons <- c(reasons, "num_group")
    }
    if (min_num_group < 5) { # 最小分组样本数 < 5
      reasons <- c(reasons, "min_num_group")
    }
    if (num_surv < 5) { # 生存事件数 < 5
      reasons <- c(reasons, "num_surv")
    }

    # if (marker_type == "SNP") {
    #   maf <- calculate_maf(survival_data$value)
    #   if (maf <= 0.05 || is.na(maf)) {
    #     reasons <- c(reasons, "maf")
    #   }
    # }
  } else {
    reasons <- c(reasons, "num_sample", "num_group", "min_num_group")
  }

  if (marker_type %in% c("gene", "mRNA", "lncRNA", "miRNA")) {
    mean_exp <- mean(survival_data$value, na.rm = TRUE)
    if (mean_exp < log2(0.1+1) || is.na(mean_exp)) {
      reasons <- c(reasons, "mean_exp")
    }
  }
  
  ## 表达阈值改为中位表达 < 0.01，以去除表达极低的RNA
  if (marker_type %in% c("miRNA")) {
    med_exp <- median(survival_data$value, na.rm = TRUE)
    if (2^med_exp - 1 < 0.01 || is.na(med_exp)) {
      reasons <- c(reasons, "med_exp")
    }
  } else if (marker_type %in% c("gene", "mRNA", "lncRNA")) {
    med_exp <- median(survival_data$value, na.rm = TRUE)
    if (med_exp < 0.01 || is.na(med_exp)) {
      reasons <- c(reasons, "med_exp")
    }
  }
  
  ## APA 过滤标准差
  if (marker_type == "APA") {
    std <- sd(survival_data$value, na.rm = TRUE)
    if (std < 0.05 || is.na(std)) {
      reasons <- c(reasons, "std")
    }
  }
  
  # med <- median(survival_data$value, na.rm = TRUE)
  # if (med < 0.1 || is.na(med)) {
  #   reasons <- c(reasons, "med")
  # }
  
  ## 去除性染色体的AS事件
  if (marker_type == "AS") {
    gene_name <- sub("_.*", "", marker)
    sex_gene_list <- scan("/home/wuzj/survival/01.data/AS/gene_in_chrXY",
                          what = "characater", quiet = TRUE)
    if (gene_name %in% sex_gene_list) {
      reasons <- c(reasons, "sex_chr")
    }
  }
  
  ## 甲基化位点过滤
  if (marker_type == "methylation") {
    probe <- marker
    ## 去除性染色体的探针
    sex_chr_probe <- scan("/home/wuzj/survival/01.data/methylation/probeMap/probe_in_chrXY",
                          what = "characater", quiet = TRUE)
    ## 去除交叉反应的探针
    cross_reactive_probe <- scan("/home/wuzj/survival/01.data/methylation/probeMap/cross_reactive_probe",
                                 what = "characater", quiet = TRUE)
    ## 去除多态性的探针
    polymorphic_probe <- scan("/home/wuzj/survival/01.data/methylation/probeMap/polymorphic_probe",
                              what = "characater", quiet = TRUE)
    if (probe %in% sex_chr_probe) {
      reasons <- c(reasons, "sex_chr")
    }
    
    if (probe %in% cross_reactive_probe) {
      reasons <- c(reasons, "cross_reactive")
    }
    
    if (probe %in% polymorphic_probe) {
      reasons <- c(reasons, "polymorphic")
    }
  }
  
  if (length(reasons) > 0) {
    # print(paste(marker, "filtered by", paste(reasons, collapse = ", ")))
    return(list(pass = FALSE, reason = reasons))
  } else {
    # print(paste(marker, "passed"))
    return(list(pass = TRUE, reason = NA))
  }
}

##### 测试 #####
# 标志物类型：gene(mRNA,lncRNA),miRNA,methylation,AS,APA,SNP,CNV,protein
# cancer_type <- "CHOL"
# marker_type <- "methylation"
# survival_type <- "DFI"
# fpath.marker <- "/home/wuzj/tmp/kmplot_data/methylation/CHOL.txt"
# fpath.clinic <- "/home/wuzj/survival/01.data/clinic/CHOL.clinic.tsv"
# output_dir <- "/home/wuzj/survival/01.data/READ/filter_res_step1"

#### 读入数据 ####
## 表型数据（有header）
header <- scan(paste0("/home/wuzj/tmp/kmplot_data/", marker_type, "/header.", cancer_type),
               what = "character", quiet = T, sep = "\t")

data.marker <- fread(fpath.marker, header = F, sep = "\t", data.table = F) %>%
  set_names(header)
  # %>% filter(ID == "rs11751451:168413393:C:T:-")
  

data.marker <- process_marker_data(data.marker, marker_type)

## 临床数据
data.clinic <- fread(fpath.clinic, header = T, sep = "\t")

##### 运行 #####
## 初始化
# 过滤情况数量统计
counts <- list(
  num_marker = length(data.marker$marker),
  passed = 0,
  missing_rate = 0,
  num_sample = 0,
  num_group = 0,
  min_num_group = 0,
  num_surv = 0,
  # maf = 0,
  med_exp = 0,
  # med = 0,
  std = 0,
  sex_chr = 0,
  cross_reactive = 0,
  polymorphic = 0
)

# 过滤后的marker列表
results <- data.frame(marker = character())

## 循环处理
for (i in 1:nrow(data.marker)) {
  
  # res <- filter_data(marker_type, data.marker[i,], data.clinic)
  res <- tryCatch({
    filter_data(marker_type, data.marker[i,], data.clinic)
  }, error = function(e) {
    ## 排错：提示并跳过
    message(paste("Error:", data.marker[i,]$marker))
    message(e)
    return(NULL)
  })
  
  if (is.null(res)) {
    next
  }
  
  if (!res$pass) {
    for (reason in res$reason) {
      counts[[reason]] <- counts[[reason]] + 1
    }
  } else {
    counts$passed <- counts$passed + 1
    results <- rbind(results, data.frame(marker = data.marker[i, "marker"]))
  }
}

# 打印结果统计
# print(paste("Number of passed markers:", counts$passed))
# for (reason in names(counts)) {
#   if (reason != "passed" && counts[[reason]] > 0) {
#     print(paste("Number of markers filtered by", reason, ":", counts[[reason]]))
#   }
# }

filter_stats <- as.data.frame(t(unlist(counts)), stringsAsFactors = FALSE)

# 保存
dir.create(output_dir, recursive = T, showWarnings = F)

fwrite(results, file.path(output_dir, paste0(basename(fpath.marker), "_filtered_marker_list")),
       sep = "\t", quote = F, row.names = F, col.names = F, na = "NA")

output_dir2 <- paste0("/home/wuzj/survival/01.data/", marker_type, "/filter_res_", survival_type)
# output_dir2 <- paste0("/home/wuzj/survival/01.data/", marker_type, "/filter_res_step1")
output_file2 <- file.path(output_dir2, paste0(cancer_type, "_filter_stats"))

if (!file.exists(output_file2)) {
  fwrite(filter_stats, output_file2, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA", append = TRUE)
} else {
  fwrite(filter_stats, output_file2, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "NA", append = TRUE)
}






