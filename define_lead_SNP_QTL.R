# September 23,2022
# Version 0
# Author: Dailu Guan
### load required packages
load_packages <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) return(TRUE)
  }else{
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}
}

### parse parameters
load_packages("argparser")
arguments <- arg_parser("*************** The script is used for finding leading SNP and defining QTL region from GWAS results ***************")
arguments <- add_argument(arguments, "--gwas", short="-g", help="path to GWAS result file, currently only support GCTA format")
arguments <- add_argument(arguments, "--is_fdr", short="-f", default='TRUE', help="Wheter FDR is done, and a clomun name: fdr, included in header")
arguments <- add_argument(arguments, "--is_signif", short="-s", default='TRUE', help="Wheter only using significant SNPs for defining QTL region")
arguments <- add_argument(arguments, "--chr", short="-c", default='Chr', help="Chromosome column")
arguments <- add_argument(arguments, "--bp", short="-b", default='bp', help="SNP position column")
arguments <- add_argument(arguments, "--pval", short="-p", default='p', help="GWAS P value column")
arguments <- add_argument(arguments, "--out_prefix", short="-o",default='QTL_region', help="GWAS P value column")
argv <- parse_args(arguments)

### Options
Chr=argv$Chr
bp=argv$bp
p=argv$pval

### Read data
load_packages("data.table")
load_packages("tidyverse")
df = fread(argv$gwas, header=T)
if (argv$is_fdr == "FALSE"){
  df$fdr = p.adjust(df[,p], method="fdr")
  fdr="fdr"
}else{
  fdr="fdr"
}
df$log10 = -log10(df[,p])

### Chromosomes with QTLs
signif_chrs=as.vector(unique(subset(df, df[,fdr] < 0.05)[,Chr]))

#### define lead snps
lead_snp_df =data.frame()
for (c in signif_chrs){
  signif_df = subset(df, df[,Chr] == c & df[, fdr] < 0.05 )
  lead_snp_df = rbind(lead_snp_df, signif_df[which.min(signif_df[,p])])
  while (nrow(signif_df) > 0){
    signif_df[which.min(signif_df[, p])] -> lead_snp_df_tmp
        left_lead = lead_snp_df_tmp[, bp] - 1000000
    if (left_lead < 0){
            left_lead = 0
    }else{
        left_lead=left_lead
    }
    right_lead = lead_snp_df_tmp[, bp] + 1000000
    signif_df = subset(signif_df, !(signif_df[,bp] >= left_lead & signif_df[,bp] <= right_lead))
    signif_df[which.min(signif_df[,p])] -> lead_snp_df_tmp
    lead_snp_df = rbind(lead_snp_df, lead_snp_df_tmp)
  }
}

res_df = data.frame()
#### define QTL regions
if (nrow(lead_snp_df) > 0){
  for (n in 1:nrow(lead_snp_df)){
    chr = as.integer(lead_snp_df[n,Chr])
    lead_pos = as.integer(lead_snp_df[n,bp])
    left_lead = as.integer(lead_snp_df[n,bp]) - 1000000
    if (left_lead < 0){
        left_lead = 0
    }else{
        left_lead=left_lead
    }
    right_lead = as.integer(lead_snp_df[n,bp]) + 1000000
    if (argv$is_signif == "TRUE"){
        left_df = subset(df, df[,Chr] == chr & df[,bp] >= left_lead & df[,bp] <= as.integer(lead_snp_df[n, bp]) & df[,fdr] <= 0.05) %>% arrange(desc({{bp}}))
        right_df = subset(df, df[,Chr] == chr & df[,bp] <= right_lead & df[,bp] >= as.integer(lead_snp_df[n, bp]) & df[,fdr] <= 0.05) %>% arrange(desc({{bp}}))
    }else{
        left_df = subset(df, df[,Chr] == chr & df[,bp] >= left_lead & df[,bp] <= as.integer(lead_snp_df[n, bp])) %>% arrange(desc({{bp}}))
        right_df = subset(df, df[,Chr] == chr & df[,bp] <= right_lead & df[,bp] >= as.integer(lead_snp_df[n, bp])) %>% arrange(desc({{bp}}))
    }
    top_pval = as.numeric(lead_snp_df[n,"log10"])

    qtl_left=as.integer(min(left_df[top_pval - left_df$log10 < 4, bp]))
    left_dis = lead_pos - qtl_left
    if (left_dis <= 500000){
        qtl_left=qtl_left
    }else{
        qtl_left=subset(left_df, left_df[,bp] >= lead_pos - 500000) %>% pull({{bp}}) %>% min 
    }
    qtl_right=as.integer(max(right_df[top_pval - right_df$log10 < 4, bp]))
    right_dis = qtl_right - lead_pos
    if (left_dis <= 500000){
        qtl_right=qtl_right
    }else{
        qtl_right=subset(right_df, right_df[,bp] <= lead_pos + 500000) %>% pull({{bp}}) %>% max
    }
    res_df = rbind(res_df, cbind(lead_snp_df[n,], qtl_left, qtl_right))
  }
}else{
  cat("No lead SNPs found...")
}
#### output
output_file=paste0(argv$out_prefix,".txt")
fwrite(res_df, output_file, sep = "\t", row.names=F, quote=F)
