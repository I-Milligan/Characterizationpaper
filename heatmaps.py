from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# 1. Load alignment
# -----------------------------
alignment_file = "y:\Documents\Python\Iris PHD\Plots and graphs\mmp2_muscle.fasta"   # change for MMP9
alignment = AlignIO.read(alignment_file, "fasta")
# -----------------------------
# 2. Calculate pairwise distances
# Poisson model approximated using identity -> distance
# -----------------------------
calculator = DistanceCalculator("blosum62")
dm = calculator.get_distance(alignment)

# Convert to dataframe
distance_matrix = pd.DataFrame(
    np.array(dm),
    index=dm.names,
    columns=dm.names
)

#print(distance_matrix)
#distance_matrix.to_csv("/Users/richardmilligan/Documents/Python/Iris PHD/Plots and graphs/MMP2_distance_matrix.csv") 
# -----------------------------
# 3. Heatmap
# -----------------------------
plt.figure(figsize=(10,8))
sns.heatmap(
    distance_matrix,
    cmap="coolwarm",
    square=True,
    linewidths=0.5,
    cbar_kws={"label": "Amino acid substitutions per site"}
)
plt.title("Pairwise Amino Acid Divergence (MMP2)")
plt.tight_layout()
plt.savefig("y:\Documents\Python\Iris PHD\Plots and graphs\MMP2_heatmap.png", dpi=300)
plt.show()

# -----------------------------
# 4. Build Neighbor-Joining Tree
# -----------------------------
constructor = DistanceTreeConstructor()
nj_tree = constructor.nj(dm)

# Save tree
Phylo.write(nj_tree, "y:\Documents\Python\Iris PHD\Plots and graphs\MMP2_tree.nwk", "newick")

# -----------------------------
# 5. Plot Tree
# -----------------------------
plt.figure(figsize=(10,8))
Phylo.draw(nj_tree, do_show=False)
plt.title("Neighbor-Joining Tree of MMP2")
plt.savefig("y:\Documents\Python\Iris PHD\Plots and graphs\MMP2_tree.png", dpi=300)
plt.show()