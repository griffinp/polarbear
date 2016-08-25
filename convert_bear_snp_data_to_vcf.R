# This is a script to convert the custom SNP genotype format for the
# files from Liu et al. (2014) 'Population genomics reveal
# recent speciation and rapid evolutionary adaptation in 
# polar bears" doi:10.1016/j.cell.2014.03.054
# to VCF format.
# The files were downloaded from http://gigadb.org/dataset/view/id/100008/

# files look like this:

# chromo	position	ref	anc	major	minor	#major	#minor	knownEM pK-EM	EG01	WG01	EG02	EG03	EG04	WG02	EG05	WG03	WG04	WG05	EG06	WG06	WG07	WG08	WG09	WG10	WG11	WG12
# scaffold79	945	A	A	A	G	6	30	0.16672GG	GG	AG	GG	AG	GG	GG	GG	GG	GG	AA	GG	GG	GG	AG	GG	GG	AG

# cols 1, 2, 5, 6 give first 4 cols of vcf format

# count total columns

##### FUNCTIONS #####

convert_genotype <-function(genotype_string_column, maj_min){
  two_alleles <- strsplit(genotype_string_column, split='')
  first_base <- sapply(two_alleles, "[[", 1)
  #print(first_base)
  second_base <- sapply(two_alleles, "[[", 2)
  call_matrix <- matrix(0, ncol=2, nrow=length(two_alleles))
  first_call <- as.numeric(first_base==maj_min[,2])
  second_call <- as.numeric(second_base==maj_min[,2])
  full_call <- paste(first_call, "/", second_call, sep="")
  return(full_call)
}

#####################

setwd(dir = "~/Documents/Teaching/polar_bear_data")
desired_scaffold <- "scaffold37"
polar_data <- read.csv("polar_bear.pooled.snp.txt", sep="\t", 
                       stringsAsFactors=FALSE, header=TRUE)

# SUBSET by scaffold immediately
polar_data <- subset(polar_data, polar_data[,1]==desired_scaffold)

polar_col_no <- ncol(polar_data)
polar_maj_min <- as.matrix(polar_data[,5:6])
polar_ind_data <- as.matrix(polar_data[,10:polar_col_no])

colnames(polar_ind_data)
# getting just 3 individuals per pop
polar_ind_data <- polar_ind_data[,c(1,3,4,2,6,8)]

### import and process brown bear data

brown_data <- read.csv("brown_bear.pooled.snp.txt", sep="\t", 
                       stringsAsFactors=FALSE, header=TRUE)

# SUBSET by scaffold immediately
brown_data <- subset(brown_data, brown_data[,1]==desired_scaffold)


brown_col_no <- ncol(brown_data)
brown_maj_min <- as.matrix(brown_data[,5:6])
brown_ind_data <- as.matrix(brown_data[,10:brown_col_no])

colnames(brown_ind_data)
# getting just 3 individuals per pop
brown_ind_data <- brown_ind_data[,c(1,2,3,8,9,10)]

### MERGING THE POLAR BEAR AND BROWN BEAR GENOTYPE DATA

brown_to_keep <- brown_data[which(brown_data$position%in%polar_data$position),]
polar_to_keep <- polar_data[which(polar_data$position%in%brown_to_keep$position),]

brown_and_polar <- data.frame(brown_to_keep[,c(1,2,5,6,10:12,17:19)],polar_to_keep[,c(10,12,13,11,15,17)])

brown_and_polar_converted <- apply(brown_and_polar[,5:ncol(brown_and_polar)], MARGIN=2, FUN=convert_genotype,
                                   maj_min=as.matrix(brown_and_polar[,3:4]))
colnames(brown_and_polar_converted) <- paste(rep(c("brown", "polar"), each=6), colnames(brown_and_polar_converted), sep="_")

brown_and_polar_output_vcf <- cbind(brown_and_polar[,c(1,2)], c("."), brown_and_polar[,c(3,4)], c("."),
                    FILTER="PASS", c("."), FORMAT="GT",
                    brown_and_polar_converted)

colnames(brown_and_polar_output_vcf) <- c("#CHROM", "POS", "ID", "REF", "ALT",
                          "QUAL", "FILTER", "INFO", "FORMAT", colnames(brown_and_polar_converted))

preheader <- '##fileformat=VCFv4.0\n##fileDate=20160406\n##source=Liu_et_al_2014_doi:10.1016/j.cell.2014.03.054\n##phasing=unphased\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"'

write.table(preheader, file="brown_and_polar_bear_scaffold37_snps.vcf", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(brown_and_polar_output_vcf, file="brown_and_polar_bear_scaffold37_snps.vcf", append=TRUE, quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")



polar_converted <- apply(data.frame(polar_ind_data), MARGIN=2, FUN=convert_genotype, maj_min=polar_maj_min)

polar_output_vcf <- cbind(polar_data[,c(1,2)], c("."), polar_data[,c(5,6)], c("."),
                    FILTER="PASS", c("."), FORMAT="GT",
                    polar_converted)

colnames(polar_output_vcf) <- c("#CHROM", "POS", "ID", "REF", "ALT",
                          "QUAL", "FILTER", "INFO", "FORMAT", colnames(polar_converted))

preheader <- '##fileformat=VCFv4.0\n##fileDate=20160406\n##source=Liu_et_al_2014_doi:10.1016/j.cell.2014.03.054\n##phasing=unphased\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"'

write.table(preheader, file="polar_bear_scaffold37_snps.vcf", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(polar_output_vcf, file="polar_bear_scaffold37_snps.vcf", append=TRUE, quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")



brown_converted <- apply(data.frame(brown_ind_data), MARGIN=2, FUN=convert_genotype, maj_min=brown_maj_min)

brown_output_vcf <- cbind(brown_data[,c(1,2)], c("."), brown_data[,c(5,6)], c("."),
                    FILTER="PASS", c("."), FORMAT="GT",
                    brown_converted)

colnames(brown_output_vcf) <- c("#CHROM", "POS", "ID", "REF", "ALT",
                          "QUAL", "FILTER", "INFO", "FORMAT", colnames(brown_converted))

preheader <- '##fileformat=VCFv4.0\n##fileDate=20160406\n##source=Liu_et_al_2014_doi:10.1016/j.cell.2014.03.054\n##phasing=unphased\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"'

write.table(preheader, file="brown_bear_scaffold37_snps.vcf", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(brown_output_vcf, file="brown_bear_scaffold37_snps.vcf", append=TRUE, quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")




