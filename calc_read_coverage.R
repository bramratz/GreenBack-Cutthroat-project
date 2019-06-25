## Get allele depth (run outside of R on Teton)
## Run time is short enough that it's OK to do it outside of the queue submission system

#module load swset gcc vcftools
#vcftools --vcf  variants_rbtassem_0.5_maf0.05.recode.vcf --extract-FORMAT-info AD --out variants_rbtassem_0.5_maf0.05.recode

## Simplify IDs using sed
## Remove the damn commas between alleles

#sed 's/\/project\/ysctrout\/emandevi\/gbcutthroat\/bwa_assem_rbt\/aln_//g' variants_miss0.5_maf0.05_ind0.95.recode.AD.FORMAT | sed 's/\.sorted\.bam//g' | sed 's/\t/,/g' > variants_miss0.5_maf0.05_ind0.95alleledepth_simplifiedIDs.txt

## Read allele depth file into R

ad <- read.csv("variants_miss0.5_maf0.05_ind0.95alleledepth_simplifiedIDs.txt", skip=1, header=F)

perlocus <- apply(ad[,3:dim(ad)[2]], 1, sum)

a1 <- ad[,seq(3,(dim(ad)[2]-1),2)]
a2 <- ad[,seq(4,dim(ad)[2],2)]

nind <- (dim(ad)[2]/2) - 1
nloci <- dim(ad)[1]

readcount <- a1+a2
perind <- apply(readcount, 2, sum) ## Think more about how to exclude low-coverage inds for other datasets using this metric

readsperindperlocus <- mean(pply(readcount, 2, sum)/nloci)

print(readsperindperlocus)
