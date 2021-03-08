library(magrittr)

## READ DATA (a table created using GATK VariantsToTable tool)
tab <- read.table("input.tsv", h=T, sep="\t")

## SAMPLE PAIRS
mlls_combn <- grep("GT",colnames(tab),val=T) %>% combn(2) %>% t

## PREPARE OUTPUT TABLE
mlls_dist <- data.frame(mlls_combn,matrix(NA,nrow(mlls_combn),3))
names(mlls_dist) <- c("anim1","anim2","total","mismatch","intraMLL")

## REMOVE SUSPICIOUS POSITIONS (homozygous in F1 hybrids, only makes sense for table of diagnostic positions)
tab <- tab[-grep("A/A|C/C|G/G|T/T",tab$csc069_GT),]
tab <- tab[-grep("A/A|C/C|G/G|T/T",tab$csc071_GT),]

## COUNT MISMATCHES FOR SAMPLE PAIRS
for(i in seq(nrow(mlls_dist))){
  x <- tab[which(is.na(tab[,grep(mlls_dist[i,1],colnames(tab))])==F & is.na(tab[,grep(mlls_dist[i,2],colnames(tab))])==F),c(grep(mlls_dist[i,1],colnames(tab)),grep(mlls_dist[i,2],colnames(tab)))]
  print(colnames(x))
  mlls_dist[i,3] <- nrow(x)
  mlls_dist[i,4] <- which(x[,1]!=x[,2]) %>% length()
}

## EXPORT TO CSV
write.csv2(mlls_dist,"pairwise_mismatches.csv")

## BOXPLOT
boxplot(list(interMLL=(mlls_dist$mismatch/mlls_dist$total)[which(mlls_dist$intraMLL==F)],intraMLL=(mlls_dist$mismatch/mlls_dist$total)[which(mlls_dist$intraMLL==T)]))

## WILCOXON RANK SUM TEST
wilcox.test(x=(mlls_dist$mismatch/mlls_dist$total)[which(mlls_dist$intraMLL==F)],y=(mlls_dist$mismatch/mlls_dist$total)[which(mlls_dist$intraMLL==T)])
