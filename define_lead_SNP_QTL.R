# September 23,2022
# Version 0
### load required packages
load_packages <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) return(TRUE)

  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}

### parse parameters
load_packages("argparser")
arguments <- arg_parser("*************** The script is used for finding leading SNP and defining QTL region from GWAS results ***************")
arguments <- add_argument(arguments, "--gwas", short="-g", help="path to GWAS result file, currently only support GCTA format")
argv <- parse_args(arguments)

### Read data
load_packages("data.table")
load_packages("tidyverse")
df = fread(argv$gwas, header=T)
df$log10 = -log10(df$p)

### Chromosomes with QTLs
signif_chrs=as.vector(unique(subset(df, fdr < 0.05)$Chr))

#### define lead snps
lead_snp_df =data.frame()
for (c in signif_chrs){
  signif_df = subset(df, Chr == c & fdr < 0.05 )
  lead_snp_df = rbind(lead_snp_df, signif_df[which.min(signif_df$p)])
  while (nrow(signif_df) > 0){
    signif_df[which.min(signif_df$p)] -> lead_snp_df_tmp
        left_lead = lead_snp_df_tmp$bp - 1000000
    if (left_lead < 0){
            left_lead = 0
    }else{
        left_lead=left_lead
    }
    right_lead = lead_snp_df_tmp$bp + 1000000
    signif_df = subset(signif_df, !(bp >= left_lead & bp <= right_lead))
    signif_df[which.min(signif_df$p)] -> lead_snp_df_tmp
    lead_snp_df = rbind(lead_snp_df, lead_snp_df_tmp)
  }
}


#### define QTL regions
res_df = data.frame()
for (n in 1:nrow(lead_snp_df)){
    chr = as.integer(lead_snp_df[n,"Chr"])
        left_lead = as.integer(lead_snp_df[n,"bp"]) - 1000000
    if (left_lead < 0){
            left_lead = 0
    }else{
        left_lead=left_lead
    }
    right_lead = as.integer(lead_snp_df[n,"bp"]) + 1000000
    left_df = subset(df, Chr == chr & bp >= left_lead & bp <= as.integer(lead_snp_df[n,"bp"])) %>% arrange(desc(bp))
    right_df = subset(df, Chr == chr & bp <= right_lead & bp >= as.integer(lead_snp_df[n,"bp"])) %>% arrange(desc(bp))
    top_pval = as.numeric(lead_snp_df[n,"log10"])

    qtl_left=as.integer(min(left_df[top_pval - left_df$log10 < 4, "bp"]))
    qtl_right=as.integer(max(right_df[top_pval - right_df$log10 < 4, "bp"]))
    res_df = rbind(res_df, cbind(lead_snp_df[n,], qtl_left,qtl_right ))
}

#### output
output_file=paste0(str_replace(argv$gwas, ".assoc.mlma.fdr.txt",""),".qtl_resions.txt")
fwrite(res_df, output_file, sep = "\t", row.names=F, quote=F)
