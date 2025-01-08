library(dplyr)
library(ggplot2)

# Check unique values of SS column
print(unique(vcf_full$SS))

# Filter variants in the normal tissue based on FILTER and SS (check if SS == "0" is correct)
normal_variants <- vcf_full %>%
  filter(FILTER == "PASS")  # Consider removing SS filter for now if necessary

# If no rows are returned, inspect the structure of your VCF file
if (nrow(normal_variants) == 0) {
  print("No variants found after filtering. Check filtering criteria.")
} else {
  # Ensure DP is numeric
  normal_variants$DP <- as.numeric(normal_variants$DP)
  
  # Filter for variants with sufficient depth (adjust threshold if necessary)
  high_quality_normal_variants <- normal_variants %>%
    filter(DP > 5)  # Try a lower threshold if necessary
  
  # If high-quality variants are found, proceed with calculations
  if (nrow(high_quality_normal_variants) > 0) {
    # Calculate mutation density
    callable_bases <- 30e6  # Adjust this value if necessary
    mutation_density_normal <- nrow(high_quality_normal_variants) / (callable_bases / 1e6)
    print(paste("Background Mutation Density (normal):", mutation_density_normal, "mutations/Mb"))
    
    # Calculate the median depth for the background mutation level
    median_depth_normal <- median(high_quality_normal_variants$DP, na.rm = TRUE)
    print(paste("Median Depth for Normal Tissue:", median_depth_normal))
    
    # Calculate reads per million (RPM) for background mutation detection
    total_depth_normal <- sum(high_quality_normal_variants$DP, na.rm = TRUE)
    reads_per_million <- (median_depth_normal / total_depth_normal) * 1e6
    print(paste("Reads per Million (RPM) to confidently call a mutation:", reads_per_million))
    
    # Plot depth distribution
    ggplot(high_quality_normal_variants, aes(x = DP)) +
      geom_histogram(binwidth = 10, fill = "skyblue", color = "black", alpha = 0.7) +
      labs(title = "Distribution of Depth in Normal Tissue", x = "Depth", y = "Frequency")
  } else {
    print("No high-quality variants found with sufficient depth.")
  }
}
