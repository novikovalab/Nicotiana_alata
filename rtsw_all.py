import seaborn as sns
from pymsaviz import MsaViz
from collections import Counter
from Bio import AlignIO
import matplotlib.pyplot as plt

# Load your MSA file using Biopython
msa_file = "../RTSW_like.fas"
alignment = AlignIO.read(msa_file, "fasta")

# Flatten the sequences for the selected region and count unique bases
#bases = [str(record.seq)[start:end] for record in alignment]  # Extract sequences in the range
bases = [str(record.seq) for record in alignment]  # Extract sequences in the range
flattened_bases = "".join(bases)  # Flatten to a single string
base_counts = Counter(flattened_bases)  # Count the unique bases

# Remove gaps ('-') from the base count and prevent coloring
if '-' in base_counts:
    del base_counts['-']

# Output the base counts (this tells you how many different bases there are)
print("Base counts:", base_counts)

# Use seaborn's color palette for more attractive colors
color_palette = sns.color_palette("husl", len(base_counts))  # Use 'husl' palette with number of unique bases

# Assign colors to each base, making sure that '-' is not colored
color_scheme = {base: color for base, color in zip(base_counts.keys(), color_palette)}

# Create MSA visualization using pymsaviz
mv = MsaViz(msa_file, show_grid=True, show_consensus=True, wrap_length=150)  # 1-based index for pymsaviz

# Set the custom color scheme for the bases (no color for '-')
mv.set_custom_color_scheme(color_scheme)

# Plot the MSA visualization
fig = mv.plotfig()

# Save the figure as a PDF
fig.savefig("rtsw_150aa.pdf", format="pdf")
