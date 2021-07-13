from activity_assay import make_path
from os import lseek
import pandas as pd
import matplotlib.pyplot as plt
import csv
import importlib
importlib.import_module("activity_assay")

# plot and save graphs
def plot():
    # print(df_amp.head())
    make_path(path, "Plots")
    for n in col:
        df1 = df_amp.loc[20:95, headers[n - 1]]
        df2 = df_deriv.loc[20:95, headers[n - 1]]
        
        plt.figure();
        df1.plot();
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='xx-small')
        plt.savefig(f"{path}Plots/{n}_{prot_substrate(n)}_t{str(trial(n))}", 
                    bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
        plt.xlabel('Temperature')
        plt.close();

        plt.figure();
        df2.plot();
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='xx-small')
        plt.savefig(f"{path}Plots/{n}_{prot_substrate(n)}_t{str(trial(n))}_derivs", 
                    bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
        plt.xlabel('Temperature')
        plt.close()
    # plt.show();

def melting():
    with open(path + 'tM_values.csv', 'w') as file:
        writer = csv.writer(file)
        for n in col:
            names = []
            row = []
            for name in headers[n - 1]:
                s = df_deriv.loc[20:95, name]
                names.append(f"{name}_{prot_substrate(n)}_t{trial(n)}")
                row.append(s.idxmin())
                # print(f"{name}_{molarity(n)}uM_t{trial(n)}: {s.idxmin()}")
            writer.writerow(names)
            writer.writerow(row)

# reading the data into panda DataFrames
path = "/Users/carlho/Documents/Shakh Lab/2021_07_12_ligand_binding_TSA/"
data = lambda fname : pd.read_excel(path + fname, 
                      sheet_name="FRET", 
                      header=0,
                      index_col=0,
                      usecols="B:NV")
df_amp = data("badkar_2021-07-12 16-38-25_CT010451 -  " + 
                "Melt Curve Amplification Results.xlsx")
df_deriv = data("badkar_2021-07-12 16-38-25_CT010451 -  " + 
                "Melt Curve Derivative Results.xlsx")

# wanted header names in DataFrames
rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
cols = range(1, 25)
headers = []

for num in cols:
    headers.append([])
    for letter in rows:
        headers[-1].append(letter + str(num))
# print(headers)
col = [n for n in range(1, 17)]  # columns where the samples are stored
# molarity = lambda n : "2" if n < 3 else "4"
# trial = lambda n : n if n < 3 else n - 3

prot_key = {
    1 : "117G",
    2 : "164G", 
    3 : "117-164G",
    4 : "WT", 
    0 : "WT"
}
trial = lambda n : ((n - 1) % 2) + 1
prot_substrate = lambda n : prot_key[int(((n % 8) + 1) / 2)] + ("_AMP" if n < 9 else "_ATP")

# print(df_amp.index)

# plot();  # uncomment to plot graphs
melting();