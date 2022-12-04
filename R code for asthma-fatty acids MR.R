


library(TwoSampleMR)
library(readxl)
library(readr)
library(MRPRESSO)
library(xlsx)
library(tidyverse)
library(readxl)
library(tidyverse)
library(ggsci)
library(UniprotR)



ao <- available_outcomes() 



met_project <-  ao[grepl("met-d", ao$id), ]


{#initial analysis 这个位置只需要改文件集的名字，别的都不用改
  
  result <- tibble(no = met_project$id, trait = met_project$trait, res = 2,
                   nsps = NA, 	b = NA, p = NA, se = NA,  r = NA, n = NA, F = NA, 
                   sz = met_project$sample_size,OR = NA, OR_L = NA, OR_U = NA, B_L =NA, B_U = NA) #
  
  n <- 0
  
}



for (i in  result$no[which(result$res == 2)]){#meta_project$id result$no[which(is.na(result$se) == TRUE)]) { #c(meta_project$id[which(meta_project$trait %in% candidate)],"prot-a-1494") Ahola$id
  exposure_dat <- extract_instruments( i,
                                       p1 = 5*10^-8,
                                       clump = FALSE,
                                       p2 = 5*10^-8)#提取暴露的summary，如果是QTL的话clump要选FALSE，并且不要设置P值的cutoff；kb是相邻SNP之间的间隔距离
  
 
  
  if(length(is.na(exposure_dat)) != 0){#为了避免exposure_dat没有SNP报错
    
    
    
    
    exposure_dat <- clump_data(exposure_dat, clump_r2=0.001, pop = "EUR")#clump，r2按照需要改
    
    
    
    
    
    
    
    
    
    
    
    out_dat <- extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = "ebi-a-GCST90014325",
      #"finn-b-N14_MALEINFERT","finn-b-CD2_MULTIPLE_MYELOMA_PLASMA_CELL_EXALLC",# "ieu-b-4957", #ieu-b-4878 ieu-b-4957"(MM 6\7\8), # 'finn-b-D3_OTHERAPLASTICANAEMIA',#finn-b-D3_OTHERAPLASTICANAEMIA finn-b-D3_AIHA_OTHER ebi-a-GCST90014325(asthma)
      proxies = TRUE,
      rsq = 0.8,
      maf_threshold = 0.3
    )
    
    
    
    if(is.null(out_dat) == FALSE ){
      dat <- harmonise_data(
        exposure_dat = exposure_dat, 
        outcome_dat = out_dat
      )
      
      
      dat$r <-2*(1-dat$eaf.exposure)*(dat$eaf.exposure)*(dat$beta.exposure)^2/(2*(1-dat$eaf.exposure)*(dat$eaf.exposure)*(dat$beta.exposure)^2+2*(1-dat$eaf.exposure)*(dat$eaf.exposure)*(dat$se.exposure)*(result$sz[which(result$no == i)])*(dat$beta.exposure)^2)
      
      
      
      res <- mr(dat)
      
      if(nrow(res) != 0){
        if( length(which(res$method == "Inverse variance weighted"))  == 0){
          result[which(result$no == i),4] <- res$nsnp[which(res$method == "Wald ratio")]
          result[which(result$no == i),5] <- res$b[which(res$method == "Wald ratio")]
          result[which(result$no == i),6] <- res$pval[which(res$method == "Wald ratio")]
          result[which(result$no == i),7] <- res$se[which(res$method == "Wald ratio")]
          result[which(result$no == i),8] <- sum(dat$r)
          result[which(result$no == i),9] <-  nrow(out_dat)
          result[which(result$no == i),10] <- ((result$sz[which(result$no == i)]-result$n[which(result$no == i)] - 1)/result$n[which(result$no == i)])*result$r[which(result$no == i)]/(1-result$r[which(result$no == i)])
          
          resor <- generate_odds_ratios(res)
          result[which(result$no == i),12] <- resor$or[which(resor$method == "Wald ratio")]
          
          result[which(result$no == i),13] <- resor$or_lci95[which(resor$method == "Wald ratio")]
          result[which(result$no == i),14] <- resor$or_uci95[which(resor$method == "Wald ratio")]
          
          result[which(result$no == i),15] <- resor$lo_ci[which(resor$method == "Wald ratio")]
          result[which(result$no == i),16] <- resor$up_ci[which(resor$method == "Wald ratio")]
          
          
        }else{
          result[which(result$no == i),4] <- res$nsnp[which(res$method == "Inverse variance weighted")]
          result[which(result$no == i),5] <- res$b[which(res$method == "Inverse variance weighted")]
          result[which(result$no == i),6] <- res$pval[which(res$method == "Inverse variance weighted")]
          result[which(result$no == i),7] <- res$se[which(res$method == "Inverse variance weighted")]
          result[which(result$no == i),8] <-sum(dat$r)
          result[which(result$no == i),9] <-  nrow(out_dat)
          result[which(result$no == i),10] <- ((result$sz[which(result$no == i)]-result$n[which(result$no == i)] - 1)/result$n[which(result$no == i)])*result$r[which(result$no == i)]/(1-result$r[which(result$no == i)])
          
          resor <- generate_odds_ratios(res)
          result[which(result$no == i),12] <- resor$or[which(resor$method == "Inverse variance weighted")]
          
          result[which(result$no == i),13] <- resor$or_lci95[which(resor$method == "Inverse variance weighted")]
          result[which(result$no == i),14] <- resor$or_uci95[which(resor$method == "Inverse variance weighted")]
          
          result[which(result$no == i),15] <- resor$lo_ci[which(resor$method == "Inverse variance weighted")]
          result[which(result$no == i),16] <- resor$up_ci[which(resor$method == "Inverse variance weighted")]
        }
        
        if(result$p[which(result$no == i)] < 0.05/nrow(result)){
          
          result[which(result$no == i),3] <- 1
          
        }else{
          result[which(result$no == i),3] <- 0
        }
      }else{
        result[which(result$no == i),3] <- 0
      }
      
      
      
      
      
      
      
    }else{
      result[which(result$no == i),3] <- 0
    }
  }else{
    result[which(result$no == i),3] <- 0
  }
  #result <- result[which(result$trait %in%c( candidate,"Interleukin-1 alpha")),]
  
  
  result$q <- p.adjust(result$p, method = "fdr", n = length(result$p))
  
  
  
  
  write_csv2(result, file = "intake_immune.txt") #"met_ap.txt")
  n <- n+1
  
  print(n)
}