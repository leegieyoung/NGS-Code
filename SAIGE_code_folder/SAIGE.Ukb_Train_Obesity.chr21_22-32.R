#!/bin/Rscript
TRAIT1 <- Sys.getenv("TRAIT1")
CHR <- Sys.getenv("CHR")
your_data_table <- read.table(paste0("/mnt/nas/gylee/0.GWAS/2.plink_result/",TRAIT1 ,"/chr",CHR,"/chr",CHR,".imputation_g_m_maf_hwe_bfile.bim"), head=F)
total_rows <- nrow(your_data_table)

# 64개로 데이터를 나눌 그룹 개수 계산
group_size <- ceiling(total_rows / 32)

# 데이터를 나누기 위한 리스트 생성
data_split_list <- vector("list", length = 32)

# 데이터를 나누어서 리스트에 담기
for (i in 1:32) {
  start_row <- (i - 1) * group_size + 1
  end_row <- min(i * group_size, total_rows)
  data_split_list[[i]] <- your_data_table[start_row:end_row, ]
  #print({class(data_split_list[[i]])})
  #print({head(data_split_list[[i]])})
  #print({head(data_split_list[[i]]$V2)})
  POS <- as.character(data_split_list[[i]]$V2)
  print({class(POS)})
  writeLines(POS ,paste0("/mnt/nas/gylee/0.GWAS/2.plink_result/",TRAIT1 ,"/chr",CHR,"/chr",CHR,"-",i,".list"))
}


