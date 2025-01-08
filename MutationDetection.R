library(vcfR)
library(dplyr)
library(VariantAnnotation)

# Load the VCF file
vcf_file <- "E:/Nikhill/Bioinformatics/results/files/Galaxy24-[VarScan_somatic_on_data_19_and_data_18].vcf" 

library(vcfR)
vcf <- read.vcfR(vcf_file)  # Attempt to load the VCF file

# Summary of VCF
print(vcf)
unique(vcf_full$INFO)  # Examine INFO content

# Convert VCF to a data frame
vcf_data <- as.data.frame(vcf@fix)

# Extract relevant INFO fields (e.g., REF, TLOD)
info <- vcfR::extract_info_tidy(vcf)

# Merge INFO fields with the VCF fixed data
vcf_full <- cbind(vcf_data, info)

# View the combined data
head(vcf_full)
head(vcf_full$INFO)

# Filter somatic variants based on SS, FILTER, and DP
somatic_variants <- vcf_full %>%
  filter(FILTER == "PASS", SS == "1", DP > 10)

# View filtered variants
head(somatic_variants)

# View the filtered normal variants
head(normal_variants)


num_somatic_variants <- nrow(somatic_variants)
print(paste("Number of somatic variants:", num_somatic_variants))

genome_size <- 30e6  # Replace with the callable genome size in bases
mutation_burden <- num_somatic_variants / (genome_size / 1e6)
print(paste("Mutation burden:", mutation_burden, "mutations/Mb"))

write.csv(somatic_variants, "filtered_somatic_variants.csv", row.names = FALSE)


