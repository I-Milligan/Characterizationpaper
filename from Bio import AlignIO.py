from Bio import AlignIO
from tkinter import Tk
from tkinter.filedialog import askopenfilename, asksaveasfilename
import csv

# Select alignment file
Tk().withdraw()
filepath = askopenfilename(title="Select aligned FASTA (.fas) file")

alignment = AlignIO.read(filepath, "fasta")

# Convert sequences to dictionary
seqs = {rec.id: str(rec.seq) for rec in alignment}
length = alignment.get_alignment_length()

# --- GROUPS (match your FASTA headers exactly) ---

cetaceans = [
    "Bottlenose_Dolphi",   # shortened header in your file
    "Humpback_Whal"
]

marine = [
    "Bottlenose_Dolphi",
    "Humpback_Whal",
    "Dugong",
    "Sealion",
    "Otter",
    "Polar_Bea"
]

ungulates = [
    "Bottlenose_Dolphi",
    "Humpback_Whal",
    "Camel",
    "Cow",
    "Hippopopotamus"
]

terrestrial = [
    "Human",
    "Mouse",
    "Rabbit",
    "Dog",
    "Beaver",
    "Capybara",
    "Tasmanian_Devi"
]

# --- DOMAIN MAP (approximate alignment regions) ---

def get_domain(pos):
    if pos < 30:
        return "Signal_peptide"
    elif pos < 90:
        return "Propeptide"
    elif pos < 215:
        return "Catalytic_start"
    elif pos < 350:
        return "Fibronectin_repeats"
    elif pos < 450:
        return "Catalytic_core"
    elif pos < 520:
        return "Hinge_linker"
    elif pos < 680:
        return "Hemopexin_domain"
    else:
        return "C_terminal"

rows = []

for i in range(length):

    column = {sp: seqs[sp][i] for sp in seqs if i < len(seqs[sp])}

    cet_res = [column[s] for s in cetaceans if s in column and column[s] != "-"]
    marine_res = [column[s] for s in marine if s in column and column[s] != "-"]
    ung_res = [column[s] for s in ungulates if s in column and column[s] != "-"]
    terr_res = [column[s] for s in terrestrial if s in column and column[s] != "-"]

    cet_set = set(cet_res)
    marine_set = set(marine_res)
    ung_set = set(ung_res)
    terr_set = set(terr_res)

    cetacean_flag = ""
    marine_flag = ""
    ungulate_flag = ""

    # Cetacean-specific
    if len(cet_set) == 1 and cet_res:
        if cet_set.isdisjoint(terr_set):
            cetacean_flag = list(cet_set)[0]

    # Marine mammal shared
    if len(marine_set) == 1 and marine_res:
        marine_flag = list(marine_set)[0]

    # Ungulate shared
    if len(ung_set) == 1 and ung_res:
        ungulate_flag = list(ung_set)[0]

    if cetacean_flag or marine_flag or ungulate_flag:

        rows.append([
            i,
            cetacean_flag,
            marine_flag,
            ungulate_flag,
            ",".join(terr_set),
            get_domain(i)
        ])

# --- Choose where to save the file ---
outfile = asksaveasfilename(
    title="Save output table",
    defaultextension=".csv",
    filetypes=[("CSV file", "*.csv")]
)

with open(outfile, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Alignment_Position",
        "Cetacean_shared",
        "Marine_shared",
        "Ungulate_shared",
        "Terrestrial_residues",
        "Domain"
    ])
    writer.writerows(rows)

print("Table saved as:", outfile)