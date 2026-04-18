from Bio import AlignIO
from tkinter import Tk
from tkinter.filedialog import askopenfilename

# Open file dialog
Tk().withdraw()
filepath = askopenfilename(title="Select alignment FASTA file")

# Load alignment
alignment = AlignIO.read(filepath, "fasta")

# Define groups
cetaceans = ["Bottlenose_Dolphi", "Humpback_Whal"]
marine_mammals = ["Bottlenose_Dolphi","Humpback_Whal","Dugong","Sealion","Otter","Polar_Bea"]

# Convert sequences to dictionary
seqs = {rec.id: str(rec.seq) for rec in alignment}
length = alignment.get_alignment_length()

cetacean_unique = []
marine_candidate = []

for i in range(length):

    column = {sp: seqs[sp][i] for sp in seqs}

    dolphin = column.get("Bottlenose_Dolphi")
    whale = column.get("Humpback_Whal")

    # Find sites where dolphin and whale match but others differ
    if dolphin == whale and dolphin not in ["-"]:

        others = [column[s] for s in column if s not in cetaceans]

        if all(res != dolphin for res in others if res != "-"):
            cetacean_unique.append((i, dolphin))

    # Marine mammal shared residues
    marine_res = [column[s] for s in marine_mammals if s in column]
    terrestrial_res = [column[s] for s in column if s not in marine_mammals]

    if len(set(marine_res)) == 1:
        if marine_res[0] not in terrestrial_res:
            marine_candidate.append((i, marine_res[0]))

print("\nCetacean-specific sites (Dolphin + Whale):")
for pos, aa in cetacean_unique:
    print(f"Position {pos}: {aa}")

print("\nMarine mammal candidate sites:")
for pos, aa in marine_candidate:
    print(f"Position {pos}: {aa}")