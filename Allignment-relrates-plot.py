import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import matplotlib.pyplot as plt
from Bio import AlignIO

# -------------------------------
# Tkinter setup
# -------------------------------
Tk().withdraw()

# --- Select aligned FASTA file ---
align_file = askopenfilename(
    title="Select your aligned FASTA file",
    filetypes=[("FASTA files (*.fasta, *.fas)", "*.fasta *.fas"), ("All files", "*.*")]
)
if not align_file:
    raise ValueError("No FASTA file selected!")

# --- Read alignment ---
alignment = AlignIO.read(align_file, "fasta")

# --- Find human sequence by header ---
human_seq = None
for record in alignment:
    if "Human" in record.id:
        human_seq = record.seq
        break
if human_seq is None:
    raise ValueError("Human sequence not found in alignment!")

# --- Select Excel file with relative rates ---
rate_file = askopenfilename(
    title="Select your Excel file with relative rates",
    filetypes=[("Excel files (*.xlsx)", "*.xlsx"), ("All files", "*.*")]
)
if not rate_file:
    raise ValueError("No Excel file selected!")

df = pd.read_excel(rate_file, header=3)  # adjust header row as needed

# --- Ask which gene ---
gene = input("Enter gene (MMP9 or MMP2): ").strip().upper()

# --- Human domain definitions ---
domain_definitions = {
    "MMP9": {
        "Signal_peptide": (1, 19),
        "Propeptide": (20, 109),
        "Catalytic_domain": (110, 444),
        "Fibronectin_type_II": (216, 389),
        "Hinge_region": (445, 508),
        "Hemopexin_like": (514, 704)
    },
    "MMP2": {
        "Signal_peptide": (1, 29),
        "Propeptide": (30, 109),
        "Catalytic_domain": (110, 449),
        "Fibronectin_type_II": (200, 364),
        "Hinge_region": (365, 433),
        "Hemopexin_like": (434, 660)
    }
}

if gene not in domain_definitions:
    raise ValueError("Gene must be 'MMP9' or 'MMP2'")

human_domains = domain_definitions[gene]

# --- Map human positions to alignment indices ---
alignment_domains = {}
for name, (start, end) in human_domains.items():
    count = 0
    start_index = None
    end_index = None
    for i, aa in enumerate(human_seq, 1):  # 1-indexed
        if aa != "-":
            count += 1
        if count == start and start_index is None:
            start_index = i
        if count == end:
            end_index = i
            break
    alignment_domains[name] = (start_index, end_index)

print(f"Mapped alignment positions for {gene}:")
for k, v in alignment_domains.items():
    print(f"{k}: {v}")

# --- Plot relative rates ---
# --- Make the line transparent ---

plt.figure(figsize=(15,5))
# Restore the line below if you want the line graph
# plt.plot(df['Site No.'], df['Rel. Rate'], color='black', lw=1, label='Relative Rate')
plt.plot(df['Site No.'], df['Rel. Rate'], color='black', lw=1, label='', linestyle='None', alpha=0.7)  # Adjust alpha for transparency
# Assign colors for domains
colors = ["lightblue", "lightgreen", "salmon", "orange", "violet", "yellow"]

for (name, (start, end)), color in zip(alignment_domains.items(), colors):
    plt.axvspan(start, end, alpha=0.3, color=color, label=name.replace("_", " "))

plt.xticks(range(0, len(human_seq)+1, 50))  # Adjust x-ticks as needed
plt.xlabel("Alignment Position")
plt.ylabel("Relative Evolutionary Rate")
plt.title(f"Site-Specific Evolutionary Rates of {gene}")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# --- Save figure automatically ---
output_file = f"{gene}_Relative_Rates.png"
plt.savefig(output_file, bbox_inches='tight', dpi=300)
print(f"Figure saved as {output_file}")

# Optional: plt.show() if you want interactive display
# plt.show()