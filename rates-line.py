import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import matplotlib.pyplot as plt

# Select file
Tk().withdraw()
filepath = askopenfilename(title="Select your Excel file", filetypes=[("Excel files", "*.xlsx")])

# Read Excel and skip extra rows
df = pd.read_excel(filepath, header=3)  # third row (0-indexed) is actual headers

# Check columns
print(df.columns)  # Should show: ['Site No.', 'Rel. Rate', '#1', '#2', '#3', '#4', '#5']

# Now plot relative rates
plt.figure(figsize=(15,5))
plt.plot(df['Site No.'],df['Rel. Rate'], color='black', lw=1, label='Relative Rate')

# Example domains for MMP9
domains = {
    "Signal_peptide": (1, 30),
    "Propeptide": (31, 90),
    "Catalytic_core": (91, 215),
    "Fibronectin_repeats": (216, 350),
    "Catalytic_core2": (351, 450),
    "Hinge_linker": (451, 520),
    "Hemopexin_domain": (521, 680),
    "C_terminal": (681, 754)
}

# Shade domains
for name, (start, end) in domains.items():
    plt.axvspan(start, end, alpha=0.2, label=name)

plt.xlabel("Alignment Position")
plt.ylabel("Relative Evolutionary Rate")
plt.title("Site-Specific Evolutionary Rates of MMP9")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()