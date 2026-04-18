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
# Assign Color
# -----------------------------
def assigncolor(val):
    if val == 'A': return 'background-color:#FFB3BA'
    elif val == 'B': return 'background-color:#FFDFBA'
    elif val == 'C': return 'background-color:#FFFFBA'
    elif val == 'D': return 'background-color:#BAFFC9'
    elif val == 'E': return 'background-color:#BAE1FF'
    elif val == 'F': return 'background-color:#E6E6FA'
    elif val == 'G': return 'background-color:#F4C2C2'
    elif val == 'H': return 'background-color:#FFDAB9'
    elif val == 'I': return 'background-color:#E0FFFF'
    elif val == 'J': return 'background-color:#FFFACD'
    elif val == 'K': return 'background-color:#D8BFD8'
    elif val == 'L': return 'background-color:#E0BBE4'
    elif val == 'M': return 'background-color:#FFCCCB'
    elif val == 'N': return 'background-color:#F5DEB3'
    elif val == 'O': return 'background-color:#F0FFF0'
    elif val == 'P': return 'background-color:#F5F5DC'
    elif val == 'Q': return 'background-color:#F0F8FF'
    elif val == 'R': return 'background-color:#E6E6FA'
    elif val == 'S': return 'background-color:#FFF0F5'
    elif val == 'T': return 'background-color:#FAFAD2'
    elif val == 'U': return 'background-color:#D3FFCE'
    elif val == 'V': return 'background-color:#FFDEAD'
    elif val == 'W': return 'background-color:#E0FFFF'
    elif val == 'X': return 'background-color:#F5F5F5'
    elif val == 'Y': return 'background-color:#FFE4E1'
    elif val == 'Z': return 'background-color:#FFF5EE'
    else: return ''

# -----------------------------
# Function to build visual for 
# selected Sequence Number
#------------------------------
def processSeq(seqnum, dfsent):
# -----------
# Find Sequence and 5 columns to each side
# If there are not 5 columns left or right 
# of the Sequence Number, adjust to show 10 
# columns of Data for reference
#------------
    lastCol = len(dfsent.columns)-10
    if seqnum <= 4:
          dfcolor = dfsent.iloc[:,0:10]
    elif seqnum >= lastCol:
        dfcolor = dfsent.iloc[:,lastCol:len(dfsent.columns)]
    else:
        startp = seqnum-4
        finishp = seqnum+6
        dfcolor = dfsent.iloc[:,startp:finishp]
    #title = alignment_file + "   Protene Comparison for selected Sequence: " + str(seqnum)
    styled = style_df_for_image(
                dfcolor,
                title=alignment_file + "      Protein Comparison for Selected Sequence: " + str(seqnum),
                index_width="180px",
                col_width="20px",
                seqstr = str(seqnum),
                )
    return(styled)

#------------------------------
# Define Style of Dataframe for Image
#------------------------------
def style_df_for_image(
    df: pd.DataFrame,
    title: str = "My Table Title",
    index_width: str = "160px",
    col_width: str = "110px",
    font_family: str = "Segoe UI, Arial, sans-serif",
    font_size: str = "14px",
    seqstr: str = '0'
    ) -> pd.io.formats.style.Styler:
    """
    Styles a DataFrame for clean image export:
      - Index column has its own width and left alignment.
      - Data columns are centered and bold.
      - Column headers are hidden.
      - Title (caption) is shown above the table.
    """
    # Base styler with caption and no column headers.
    styler = (
        df.style
          .map(assigncolor,subset=[str(seqstr)])
          .hide(axis="columns")  # remove column headers
          .set_caption(title)    # add title
          .set_table_styles([
              # Table layout and typography
              {"selector": "table", "props": [
                  ("border-collapse", "collapse"),
                  ("table-layout", "fixed"),         # makes width obey CSS widths
                  ("width", "auto"),
                  ("font-family", font_family),
                  ("font-size", font_size),
                  ("color", "#111"),
                  ("background-color", "white"),
              ]},
              # Caption (title) styling
              {"selector": "caption", "props": [
                  ("caption-side", "top"),
                  ("text-align", "center"),
                  ("font-weight", "700"),
                  ("font-size", "18px"),
                  ("padding", "6px 0 10px 0"),
              ]},
              # Index column cells (row header cells)
              {"selector": "th.row_heading", "props": [
                  ("width", index_width),
                  ("text-align", "left"),
                  ("font-weight", "600"),
                  ("padding", "8px 12px"),
                  ("border", "1px solid #D0D7DE"),
                  ("background-color", "#F6F8FA"),
                  ("white-space", "nowrap"),
              ]},
              # Hide index name cell & any blank corner cells (just in case)
              {"selector": "th.index_name", "props": [("display", "none")]},
              {"selector": "th.blank", "props": [("display", "none")]},

              # Data cells (non-index cells)
              {"selector": "td", "props": [
                  ("width", col_width),
                  ("text-align", "center"),
                  ("font-weight", "700"),
                  ("padding", "8px 12px"),
                  ("border", "none !important"),
                  ("white-space", "nowrap"),
              ]},
          ])
    )

    return styler

#------------------------------
# Ask for the export filename
#------------------------------
outputfile = input("Please specify the output filename to use for the output (will be generated as and excel file): ")
outputfile = str(filepath) + "\\" + outputfile + ".xlsx"

# -----------------------------
#  Request Sequence #
# Loop until Exit is sellected
# -----------------------------
ExitProgram = False
with pd.ExcelWriter(outputfile,engine='openpyxl') as writer:
    while not ExitProgram:
        ret = input("Please Select a valide Sequence Number 0-" + str(maxseq) + " adn Q or q to exit: ")
        if ret == "Q" or ret =="q":
            ExitProgram = True
        else: 
            dftemp = processSeq(int(ret),DF1)
            dftemp.to_excel(writer, sheet_name=ret, index=True)
            html = dftemp.to_html()
 #           htmlname = str(filepath) + "\\" + alignment_file + "_Sequence" + str(ret) + ".html"
            pngname = str(filepath) + "\\" + alignment_file + "_Sequence_" + str(ret) + ".png"
 #           with open(htmlname, "w", encoding="utf-8") as f:
 #               f.write(html)
            dfi.export(dftemp, pngname)
            
