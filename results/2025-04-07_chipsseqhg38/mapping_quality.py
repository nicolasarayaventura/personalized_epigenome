import os
import sys
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

bam_path = sys.argv[1]
output = sys.argv[2]

bam = pysam.AlignmentFile(bam_path, "rb")

mapq_counts = []
mapq0_count = 0

for read in bam:
    if read.is_secondary or read.is_supplementary:
        continue  # skip non-primary alignments
    mapq = read.mapping_quality
    mapq_counts.append(mapq)
    if mapq == 0:
        mapq0_count += 1

bam.close()

print("Number of primary alignments with MAPQ=0: {}".format(mapq0_count))

# Plot
plt.hist(mapq_counts, bins=150, color='skyblue', edgecolor='black')
plt.title("Mapping Quality Distribution (Primary Alignments)")
plt.xlabel("MAPQ")
plt.ylabel("Read Count")
plt.tight_layout()
plt.savefig(output)
