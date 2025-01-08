import math
from scipy.stats import norm

# Parameters
total_reads = 1000000  # Total sequencing depth (1 million reads)
num_pmp = 10  # Number of top PMPs
reads_per_pmp = total_reads / num_pmp  # Reads per PMP
error_rate = 0.01  # Background error rate (1% error)
confidence_level = 0.99  # Desired confidence level

# Calculate the Z-score for 99% confidence (one-tailed)
z_score = norm.ppf(confidence_level)

# Calculate the threshold using normal approximation
threshold_reads = reads_per_pmp * error_rate + z_score * math.sqrt(reads_per_pmp * error_rate * (1 - error_rate))

# Print the result
print(f"Threshold of reads required for each PMP to confidently call Tissue #2: {math.ceil(threshold_reads)} reads")
