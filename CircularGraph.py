from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import toytree
import toyplot.png

# Build a Newick Tree from loaded Fasta File
# Read the aligned FASTA file
aln = AlignIO.read("mmp2-mega.fas", "fasta")

# Calculate distance matrix using identity (p-distance for proteins)
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

# Build tree using Neighbor-Joining
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Save the tree to Newick format
Phylo.write(tree, "species.nwk", "newick")

print("Tree built and saved to species.nwk")

# Load the Newick tree
tree = toytree.tree("species.nwk")

# Plot the tree in a circular layout with node labels
# that will be rendered to html

canvas, axes, mark = tree.draw(
    width=600,
    height=600,
    tree_style="c",
    node_labels="name"
)

# Save the plot to a html file
toyplot.png.render(canvas, "circular_tree.png")     

print("Circular tree plot saved to circular_tree.png")
