## 进行生存分析（10年生存）并绘制生存曲线（TCGA）
## 传入癌症类型，标志物类型，标志物数据路径和临床数据路径
## 标志物类型：gene(lncRNA),miRNA,methylation,AS,APA,SNP,CNV,protein
## usage：Rscript do_surv.R cancer_type marker_type fpath.marker fpath.clinic output_dir survival_type

#### 导包 ####
suppressMessages({
  library(data.table)
  library(tidyverse)
  library(survival)
  library(survminer)
  library(ggsci)
})

#### 传入参数 ####
args<-commandArgs(TRUE)
### 癌症类型
cancer_type <- args[1]
### 标志物类型
marker_type <- args[2]
### 标志物数据路径
fpath.marker <- args[3]
### 临床数据路径
fpath.clinic <- args[4]
### 输出路径
output_dir <- args[5]
### 作图路径
# plot_dir <- args[6]
### 结局类型
survival_type <- args[6]

#### 数据处理 ####
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

filter_samples <- function(data.marker) {
  
  # 去除样本ID的多余后缀
  # colnames(data.marker) <- sub("^([^\\-]+-[^\\-]+-[^\\-]+-[^\\-]+).*", "\\1", colnames(data.marker))
  
  # 去除癌旁样本
  sample_ids <- colnames(data.marker)[grep("^TCGA", colnames(data.marker))]
  rm_sample_list <- grep("^TCGA-\\w{2}-\\w{4}-1", sample_ids, value = TRUE)
  
  # 每个患者仅保留一个分割样本的数据
  non_a_sample_list <- grep("^TCGA-\\w{2}-\\w{4}-\\d{2}[^A]", sample_ids, value = TRUE)
  
  
  for (sample_id in non_a_sample_list) {
    patient_id <- sub("-[0-9]{2}[A-Z]{1}.*$", "", sample_id)
    
    patient_samples <- grep(paste0("^", patient_id, "-"), sample_ids, value = TRUE)
    sorted_samples <- sort(patient_samples)
    rm_sample_list <- c(rm_sample_list, sorted_samples[-1])
  }
  
  # 去除复发癌症样本
  # recurrent_samples <- grep("^TCGA-\\w{2}-\\w{4}-02.*", sample_ids, value = TRUE)
  # rm_sample_list <- c(rm_sample_list, recurrent_samples)
  
  # 去除非原发癌（需改进）
  primary_sample_codes <- c("01", "03")
  non_primary_samples <- sample_ids[sapply(sample_ids, function(id) {
    parts <- strsplit(id, "-")[[1]]  # 拆分样本ID
    if (length(parts) >= 4) {
      sample_type <- substr(parts[4], 1, 2)  # 获取第四部分的前两个字符
      return(!(sample_type %in% primary_sample_codes))  # 返回是否为非原发癌样本
    }
    return(FALSE)  # 对于较短的ID直接返回FALSE
  })]
  rm_sample_list <- c(rm_sample_list, non_primary_samples)
  
  # 去除不需要的样本列
  data.marker_filtered <- dplyr::select(data.marker, -dplyr::one_of(rm_sample_list))
  return(data.marker_filtered)
}

process_survival_data <- function(marker_type, outcome_type, data.marker, data.clinic) {
  ### 生成生存分析数据
  ## 宽格式转长格式，用于与临床数据合并
  tmp <- gather(data.marker, key = sampleID)[-1,]
  
  ## sampleID 处理
  tmp <- tmp %>%
    mutate(sampleID = sub("^([^\\-]+-[^\\-]+-[^\\-]+).*", "\\1", sampleID)) %>%
    mutate(sampleID = gsub("_", "-", sampleID),
           value = as.numeric(value))
  
  ## 合并表型数据与临床数据
  merge_data <- merge(tmp, data.clinic, by="sampleID") %>%
    filter(!is.na(value))
  
  ## 分组
  if (marker_type == "SNP") {
    merge_data <- mutate(merge_data, group = factor(value, levels = c(0, 1, 2))) %>%
      mutate(group = droplevels(group))  # 移除未使用的级别
  } else if (marker_type == "CNV") {
    merge_data <- mutate(merge_data, group = factor(
      ifelse(value == 0, "normal", ifelse(value > 0, "gain", "loss")),
      levels = c("normal", "gain", "loss")
    )) %>%
      mutate(group = droplevels(group))  # 移除未使用的级别
  } else {
    merge_data <- mutate(merge_data, group = factor(
      ifelse(value > median(merge_data$value, na.rm = TRUE), "high", "low"),
      levels = c("low", "high")
    )) %>%
      mutate(group = droplevels(group))  # 移除未使用的级别
  }
  
  ## 动态选择生存结局类型
  time_col <- paste0(outcome_type, ".time")
  status_col <- outcome_type
  
  survival_data <- merge_data %>%
    mutate(status = .data[[status_col]]) %>%
    filter(!(is.na(value) | is.na(.data[[time_col]]) | is.na(.data[[status_col]]))) %>%
    mutate(status = ifelse(.data[[time_col]] > 3650, 0, status)) %>%
    mutate(time = ifelse(.data[[time_col]] > 3650, 3650, .data[[time_col]])) %>%
    mutate(time = time / 365)
  
  return(survival_data)
}

#### 生存分析 ####
do_survplot <- function(survival_data, marker, marker_type, cancer_type) {
  ### 拟合生存曲线（KM）
  fit.surv <- survfit(Surv(time, status) ~ group, data = survival_data)
  
  ## 标注基因型对应的样本数量
  sample_counts <- table(survival_data$group)
  
  ## 处理marker名称
  if (marker_type == "gene") {
    marker <- sub(".*\\|", "", marker)
  } else if (marker_type == "SNP") {
    marker <- sub(":.*", "", marker)
  }
  
  ## 标签
  if (marker_type == "SNP") {
    labels <- sapply(names(sample_counts), function(name) {
      count <- sample_counts[name]
      genotype_label <- ifelse(name == "0", "AA", ifelse(name == "1", "Aa", ifelse(name == "2", "aa", name)))
      paste(genotype_label, " (", count, ")", sep = " ")
    })
  } else {
    labels <- sapply(names(sample_counts), function(name) {
      count <- sample_counts[name]
      paste(name, "(", count, ")", sep = "")
    })
  }
  
  ## 如果是CNV，重新排序group，并调整标签顺序
  if (marker_type == "CNV") {
    new_order <- c("gain", "lstatuss", "normal")
    
    # 转换为因子，确保分组顺序一致
    survival_data$group <- factor(survival_data$group, levels = new_order)
    
    # 重新排序标签
    labels <- labels[new_order]
  }
  
  ## 调色板
  if (length(unique(survival_data$group)) == 2) {
    color_palette <- c("#E91E24", "#303094")  # 红蓝
  } else {
    color_palette <- c("#E91E24", "#303094", "#FFA500")  # 红蓝橙
  }
  
  ### 绘制生存曲线
  plt <- ggsurvplot(
    fit.surv, data = survival_data,
    pval = TRUE,
    censor.shape = "+", censor.size = 3,
    legend.title = "", 
    legend.labs = labels,
    palette = color_palette,
    xlim = c(0, 10), break.time.by = 2,
    xlab = "Time (years)", ylab = "Survival Probability",
    title = paste("KM Plot for", marker, "in", cancer_type),
    ggtheme = theme_survminer()
  )
  
  plt <- plt$plot + theme(plot.title = element_text(hjust = 0.5),
                          legend.direction = "horizontal")
  
  plt
  return(plt)
}

do_survival <- function(survival_data, marker, marker_type, cancer_type) {
  if(length(unique(survival_data$group)) > 1) {
    
    ### LogRank检验
    km_fit <- survdiff(Surv(time, status) ~ group, data = survival_data)
    # km_fit <- survdiff(Surv(OS.time, OS) ~ group, data = survival_data)
    
    ### UniCox
    unicox_model <- coxph(Surv(time, status) ~ group, data=survival_data, na.action=na.exclude)
    unicox_res <- summary(unicox_model)
    
    ### MultiCox
    ## 处理单一性别和肿瘤分期缺失的情况
    switch_case <- function(marker, cancer_type, survival_data) {
      if (!all(is.na(survival_data$stage)) && length(unique(survival_data$gender))!=1) {
        formula <- as.formula(paste("Surv(time, status) ~ group + age + stage + gender"))
      } else if (all(is.na(survival_data$stage)) && length(unique(survival_data$gender))==1) {
        cat("stage and gender is excluded in multicox for", marker, "in", cancer_type, "\n")
        formula <- as.formula(paste("Surv(time, status) ~ group + age"))
      } else if (all(is.na(survival_data$stage))) {
        cat("stage is excluded in multicox for", marker, "in", cancer_type, "\n")
        formula <- as.formula(paste("Surv(time, status) ~ group + age + gender"))
      } else if (length(unique(survival_data$gender))==1) {
        cat("gender is excluded in multicox for", marker, "in", cancer_type, "\n")
        formula <- as.formula(paste("Surv(time, status) ~ group + age + stage"))
      } 
      return(formula)
    }
    
    formula <- switch_case(marker, cancer_type, survival_data)
    multicox_model <- tryCatch({
      coxph(formula, data = survival_data, na.action = na.exclude)
    }, error = function(e) {
      cat("Error occurred in multicox for", marker, "in", cancer_type, "\n")
      NULL
    })
    
    multicox_res <- summary(multicox_model)
    
    ## 保存结果
    if (marker_type == "SNP") {
      res <- data.frame(
        marker = marker,
        
        KMp = 1 - pchisq(km_fit$chisq, df=length(levels(factor(survival_data$group)))-1),
        
        HR.unicox = unicox_res$coefficients[1, "exp(coef)"],
        L95CI.unicox = unicox_res$conf.int[1, "lower .95"],
        H95CI.unicox = unicox_res$conf.int[1, "upper .95"],
        UniCoxp = unicox_res$coefficients[1, "Pr(>|z|)"],
        
        HR.multicox = if (!is.null(multicox_model)) multicox_res$coefficients[1, "exp(coef)"] else NA,
        L95CI.multicox = if (!is.null(multicox_model)) multicox_res$conf.int[1, "lower .95"] else NA,
        H95CI.multicox = if (!is.null(multicox_model)) multicox_res$conf.int[1, "upper .95"] else NA,
        MultiCoxp = if (!is.null(multicox_model)) multicox_res$coefficients[1, "Pr(>|z|)"] else NA,
        
        N = length(unique(survival_data$sampleID)),
        G0 = median(survival_data$time[survival_data$group == "0"], na.rm = TRUE),
        G1 = median(survival_data$time[survival_data$group == "1"], na.rm = TRUE),
        G2 = median(survival_data$time[survival_data$group == "2"], na.rm = TRUE),
        AA = sum(survival_data$group == "0"),
        Aa = sum(survival_data$group == "1"),
        aa = sum(survival_data$group == "2")
      )
    } else if (marker_type == "CNV") {
      res <- data.frame(
        marker = marker,
        
        KMp = 1 - pchisq(km_fit$chisq, df=length(levels(factor(survival_data$group)))-1),
        
        HR.unicox = unicox_res$coefficients[1, "exp(coef)"],
        L95CI.unicox = unicox_res$conf.int[1, "lower .95"],
        H95CI.unicox = unicox_res$conf.int[1, "upper .95"],
        UniCoxp = unicox_res$coefficients[1, "Pr(>|z|)"],
        
        HR.multicox = if (!is.null(multicox_model)) multicox_res$coefficients[1, "exp(coef)"] else NA,
        L95CI.multicox = if (!is.null(multicox_model)) multicox_res$conf.int[1, "lower .95"] else NA,
        H95CI.multicox = if (!is.null(multicox_model)) multicox_res$conf.int[1, "upper .95"] else NA,
        MultiCoxp = if (!is.null(multicox_model)) multicox_res$coefficients[1, "Pr(>|z|)"] else NA,
        
        N = length(unique(survival_data$sampleID)),
        G0 = median(survival_data$time[survival_data$group == "normal"], na.rm = TRUE),
        G1 = median(survival_data$time[survival_data$group == "loss"], na.rm = TRUE),
        G2 = median(survival_data$time[survival_data$group == "gain"], na.rm = TRUE),
        AA = sum(survival_data$group == "normal"),
        Aa = sum(survival_data$group == "loss"),
        aa = sum(survival_data$group == "gain")
      )
    } else {
      res <- data.frame(
        marker = marker,
        KMp = 1 - pchisq(km_fit$chisq, df=length(levels(factor(survival_data$group)))-1),
        HR.unicox = unicox_res$coefficients[1, "exp(coef)"],
        L95CI.unicox = unicox_res$conf.int[1, "lower .95"],
        H95CI.unicox = unicox_res$conf.int[1, "upper .95"],
        UniCoxp = unicox_res$coefficients[1, "Pr(>|z|)"],
        HR.multicox = if (!is.null(multicox_model)) multicox_res$coefficients[1, "exp(coef)"] else NA,
        L95CI.multicox = if (!is.null(multicox_model)) multicox_res$conf.int[1, "lower .95"] else NA,
        H95CI.multicox = if (!is.null(multicox_model)) multicox_res$conf.int[1, "upper .95"] else NA,
        MultiCoxp = if (!is.null(multicox_model)) multicox_res$coefficients[1, "Pr(>|z|)"] else NA,
        N = length(unique(survival_data$sampleID)),
        G0 = median(survival_data$time[survival_data$group == "high"], na.rm = TRUE),
        G1 = median(survival_data$time[survival_data$group == "low"], na.rm = TRUE),
        AA = sum(survival_data$group == "high"),
        Aa = sum(survival_data$group == "low")
      )
    }
    
    return(res)
    
  } else {
    cat(marker, "in", cancer_type, "filtered.\n")
    return(NULL)
  }
}

#### 测试 ####
# cancer_type <- "LAML"
# marker_type <- "SNP"
# survival_type <- "OS"
# fpath.marker <- "/home/wuzj/tmp/kmplot_data/SNP/LAML.txt"
# fpath.clinic <- "/home/wuzj/survival/01.data/clinic/LAML.clinic.tsv"
# output_dir <- paste0("/home/wuzj/survival/03.result/res_",survival_type,"/TGCT")
# plot_dir <- "/home/wuzj/survival/03.result/SNP/kmplot/ACC"

#### 设定路径 ####
dir.create(output_dir, recursive = T, showWarnings = F)
# dir.create(plot_dir, recursive = T, showWarnings = F)

#### 初始化结果数据框 ####
if (marker_type %in% c("SNP","CNV")) {
  results <- data.frame(
    marker = character(),
    
    KMp = double(),
    
    HR.unicox = double(),
    L95CI.unicox = double(),
    H95CI.unicox = double(),
    UniCoxp = double(),
    
    HR.multicox = double(),
    L95CI.multicox = double(),
    H95CI.multicox = double(),
    MultiCoxp = double(),
    
    N = integer(),
    G0 = double(),
    G1 = double(),
    G2 = double(),
    AA = integer(),
    Aa = integer(),
    aa = integer(),
    stringsAsFactors = FALSE
  )
} else {
  results <- data.frame(
    marker = character(),
    
    KMp = double(),
    
    HR.unicox = double(),
    L95CI.unicox = double(),
    H95CI.unicox = double(),
    UniCoxp = double(),
    
    HR.multicox = double(),
    L95CI.multicox = double(),
    H95CI.multicox = double(),
    MultiCoxp = double(),
    
    N = integer(),
    G1 = double(),
    G2 = double(),
    N1 = integer(),
    N2 = integer(),
    
    stringsAsFactors = FALSE
  )
}

#### 读入数据 ####
## 表型数据（有header）
# header <- scan(paste0("/home/wuzj/survival/01.data/", marker_type, "/header.", cancer_type),
header <- scan(paste0("/home/wuzj/tmp/kmplot_data/", marker_type, "/header.", cancer_type),
               what = "character", quiet = T, sep = "\t")

# data.marker <- fread(fpath.marker, header = F, sep = "\t", data.table = F) %>%
#   set_names(header)

data.marker <- fread(fpath.marker, header = F, sep = "\t", data.table = F) %>%
  set_names(header)

# fpath.marker <- '/home/wuzj/CLEC3A.txt'
# data.marker <- fread(fpath.marker, header = T, sep = "\t", data.table = F)


# 蛋白换数据
# data.marker <- fread(fpath.marker, header = T, sep = "\t", data.table = F) %>%
#   rename(.,"sample" = Sample_ID)

# data.marker <- filter(data.marker, sample %in% marker_list)

data.marker <- filter_samples(data.marker)

data.marker <- process_marker_data(data.marker, marker_type)

# data.marker <- data.marker %>%
#   select(!duplicated(names(data.marker)))

## 临床数据
data.clinic <- fread(fpath.clinic, header = T, sep = "\t")
# %>% rename(sampleID=demographic.submitter_id)

#### 批量处理 ####
### 一次处理一个marker
for (i in 1:nrow(data.marker)) {
  # i <- 1
  marker <- data.marker[i,]$marker
  # print(marker)
  
  tryCatch({
    survival_data <- process_survival_data(marker_type, survival_type, data.marker[i, ], data.clinic)
    res <- do_survival(survival_data, marker, marker_type, cancer_type)
    
    if (!is.null(res)) {
      results <- rbind(results, res)
      
      # 是否画图：
      # if (res$KMp < 0.05 && res$UniCoxp < 0.05 && res$MultiCoxp < 0.05) {
        # plt <- do_survplot(survival_data, marker, marker_type, cancer_type)
        # plt
      #   ## 保存KMplot
      #   ggsave(plt, filename = paste0(plot_dir,"/",cancer_type,"-",marker,".png"),
      #        width = 6, height = 5)
      # }
    }
  }, error = function(e) {
    # cat("!!!\n")
    cat("Unkown error occurred in multicox for", marker, "in", cancer_type, "\n")
    # print(traceback())
    return()
  })
}


# colnames(results) <- c("marker", "KMp", "HR.unicox", "L95CI.unicox", "H95CI.unicox",
#                        "UniCoxp", "HR.multicox", "L95CI.multicox", "H95CI.multicox",
#                        "MultiCoxp", "N", "G1", "G2", "N1", "N2")

### 保存结果
fwrite(results, file.path(output_dir, paste0(basename(fpath.marker), "_survival_result")),
       sep = "\t", quote = F, row.names = F, col.names = F, na = "NA")

# 蛋白修改
# fwrite(results, file.path(output_dir, paste0(cancer_type, "_survival_result")),
#        sep = "\t", quote = F, row.names = F, col.names = T, na = "NA")






