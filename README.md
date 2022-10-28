# QTL_definition
Rscript for finding leading SNPs and defining QTL regions

Currently the script only supports output from GCTA tool, others will be included soon.

The script adopts description from:

Bouwman, A.C., Daetwyler, H.D., Chamberlain, A.J. et al. Meta-analysis of genome-wide association studies for cattle stature identifies common genes that regulate body size in mammals. Nat Genet 50, 362–367 (2018). https://doi.org/10.1038/s41588-018-0056-5

```
A QTL was defined as a chromosomal region where adjacent pairs of significant variants were less than 1 Mb from each other. 
Within each locus, the most significant variant was taken as the lead variant. From the lead variant within such a locus, 
a more conservative QTL locus was defined on the basis of a –log10 (P value) drop off of 4, i.e., the difference between 
the –log10 (P value) of the lead variant and variants on either side moving further until all SNPs had a difference in 
–log10 (P value) from the lead SNP of greater than 4 (if the drop in –log10 (P value) was greater than 4, then decreased 
again, the procedure continued until all further SNPs had a difference in –log10 (P value) from the lead SNP of greater 
than 4). The maximum distance considered was 0.5 Mb on either side of the lead variant.
```
