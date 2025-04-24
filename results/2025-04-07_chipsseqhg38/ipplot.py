import sys
import numpy as np
import matplotlib.pyplot as plt
import pysam

file_l1 = sys.argv[1]
file_l2 = sys.argv[2]

output = sys.arv[3]

# Function to read BAM data using pysam
def read_bam_data(bam_file):
    # Open the BAM file using pysam
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    # Fetch all reads from the entire chromosome
    positions = []
    coverage = []
    
    # Fetch reads across all chromosomes (or a specific chromosome if needed)
    for read in samfile.fetch():
        positions.append(read.reference_start)  # Position of the read
        coverage.append(1)  # Each read counts as 1 unit of coverage

    samfile.close()
    
    return np.array(positions), np.array(coverage)

# Load data from both BAM files
x_l1, y_l1 = read_bam_data(file_l1)
x_l2, y_l2 = read_bam_data(file_l2)

# Aggregate coverage by position (if multiple reads map to the same position, sum their coverage)
def aggregate_coverage(x, y):
    # Get unique positions and sum the coverage at each unique position
    unique_positions, aggregated_coverage = np.unique(x, return_counts=True)
    return unique_positions, aggregated_coverage

# Aggregate the coverage data
x_l1, y_l1 = aggregate_coverage(x_l1, y_l1)
x_l2, y_l2 = aggregate_coverage(x_l2, y_l2)

# Plot IP Strength
plt.plot(x_l1, y_l1, label='L001', color='blue')
plt.plot(x_l2, y_l2, label='L002', color='green')

# Calculate the point of rapid rise (simple derivative threshold)
threshold = 1.5  # Customize the threshold for rapid rise
diff_l1 = np.diff(y_l1)
diff_l2 = np.diff(y_l2)

# Find where the difference exceeds the threshold (indicating rapid rise)
rapid_rise_x_l1 = x_l1[np.argmax(diff_l1 > threshold)]
rapid_rise_x_l2 = x_l2[np.argmax(diff_l2 > threshold)]

# Add dashed lines at the rapid rise points
plt.axvline(x=rapid_rise_x_l1, color='red', linestyle='--', label='Rapid Rise L001')
plt.axvline(x=rapid_rise_x_l2, color='orange', linestyle='--', label='Rapid Rise L002')

# Customize plot
plt.xlabel('Position on Chromosome')
plt.ylabel('IP Strength')
plt.title('IP Strength Plot with Dashed Lines Indicating Rapid Rise')
plt.legend()

# Save the plot
output_file = 'ip_plot_with_dashed_line.png'
plt.savefig(output_file)
plt.close()

# Print output for logging
print(f"Plot saved to {output_file}")
