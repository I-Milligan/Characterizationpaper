from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import dataframe_image as dfi
from pandas.io.formats.style import Styler
import tkinter as tk
from tkinter import filedialog
from pathlib import Path

# --------------------------------
# function to select file to procrss
# --------------------------------
def select_file():
    root = tk.Tk()
    root.withdraw()  # hide the empty tkinter window

    file_path = filedialog.askopenfilename(
        title="Select a file",
        filetypes=[("All files", "*.*")]
    )
    return file_path

# -----------------------------
#  Load alignment and global variables
# -----------------------------
path = select_file()
print("Selected:", path)
# ------------------------
# Find the path to use later
# ------------------------
filepath = Path(path).parent
# ------------------------
# Get filename to use in Title and naming of other files
# ------------------------
alignment_file = Path(path).stem
# ----------------------
# Align data for conversion to dataframe
# ----------------------
alignment = AlignIO.read(path, "fasta")

# -----------------------------
#  Build Index from seqID for use in Dataframe
# -----------------------------
DFList = []
DFCol = []
#-- Build Index List for DataFrame
for record in alignment:
    DFList.append(record.id)
#-- Find Max Number of Columns
maxseq = alignment.get_alignment_length()-1
#-- Build DataFrame Column list from Column Index
x=0
while x<=maxseq:
    DFCol.append(str(x))
    x=x+1
# -----------------------------
#  Load into a DataFrame
# -----------------------------
DF1 = pd.DataFrame(np.array(alignment), index=DFList, columns=DFCol)
#  DF1.to_excel("y:\Documents\Python\Iris PHD\Plots and graphs\exporttest.xlsx")

# -----------------------------
# Select Animals to Compare against the rest
# -----------------------------
def select_animals():
    print("Select the animals to compare against the rest. Type 'done' when finished.")
    print(DFList)
    selected_animals = []
    while True:
        animal = input("Enter animal name (or 'done' to finish): ")
        if animal.lower() == 'done':
            break
        elif animal in DF1.index:
            selected_animals.append(animal)
        else:
            print(f"{animal} not found in the dataset. Please try again.")
    return selected_animals

# -----------------------------
# Main
# -----------------------------
# Start by selecting animals to compare against the rest
# ------------------------------
selected_animals = select_animals()

# -----------------------------
# Find all the sequences where the selected animals match and not equal to -
# -----------------------------
matching_sequences = []
for col in DF1.columns:
    if all(DF1.loc[animal, col] == DF1.loc[selected_animals[0], col] for animal in selected_animals) and any(DF1.loc[animal, col] != '-' for animal in selected_animals):
        matching_sequences.append(col)

# ------------------------------
# Using the matching_sequences list, find all where the other animals to not match the selected animals
# ------------------------------
non_matching_sequences = []
for col in DF1.columns:
    if col in matching_sequences:
        if any(DF1.loc[:, col] != DF1.loc[selected_animals[0], col]):
            non_matching_sequences.append(col)
# print("Non-matching sequences:", non_matching_sequences)
# -----------------------------
# Write out the non_matching sequences to a CSV file
# -----------------------------
output_file = filepath / f"{alignment_file}_sequences.csv"
DF1.loc[:, non_matching_sequences].to_csv(output_file)
print(f"Non-matching sequences saved to {output_file}")
