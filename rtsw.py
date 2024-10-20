from pymsaviz import MsaViz
from collections import Counter
from Bio import AlignIO
import random
import matplotlib.pyplot as plt

# Load your MSA file using Biopython
msa_file = "./RTSW_like.fas"
alignment = AlignIO.read(msa_file, "fasta")

# Define the start and end positions for the region of interest (1-based index in fasta)
start = 902 - 1  # Convert to 0-based index for python
end = 936

# Flatten the sequences for the selected region and count unique bases
bases = [str(record.seq)[start:end] for record in alignment]  # Extract sequences in the range
flattened_bases = "".join(bases)  # Flatten to a single string
base_counts = Counter(flattened_bases)  # Count the unique bases

# Remove gaps ('-') from the base count and prevent coloring
if '-' in base_counts:
    del base_counts['-']

# Output the base counts (this tells you how many different bases there are)
print("Base counts:", base_counts)

# Generate a unique color for each base using a visually appealing color palette
color_palette = plt.cm.Set3.colors  # Use Set3 colormap from matplotlib for distinct colors
colors = list(color_palette) * (len(base_counts) // len(color_palette) + 1)  # Repeat palette if needed

# Assign colors to each base, making sure that '-' is not colored
color_scheme = {base: color for base, color in zip(base_counts.keys(), colors)}

# Create MSA visualization using pymsaviz
mv = MsaViz(msa_file, start=902, end=936, show_consensus=True)  # 1-based index for pymsaviz

# Set the custom color scheme for the bases (no color for '-')
mv.set_custom_color_scheme(color_scheme)

# Plot the MSA visualization
fig = mv.plotfig()

# Save the figure as a PDF
fig.savefig("rtsw_helix.pdf", format="pdf")
