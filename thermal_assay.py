from os import lseek
import pandas as pd
import matplotlib.pyplot as plt
import csv

# plot and save graphs
def plot():
    # print(df_amp.head())
    for n in col:
        df1 = df_amp.loc[20:95, headers[n - 1][:7]]
        df2 = df_deriv.loc[20:95, headers[n - 1][:7]]
        
        df1.plot();
        # comment/uncomment line below to save files
        plt.savefig(f"{path}Plots/{molarity(n)}uM_t{str(trial(n))}")
        df2.plot();
        # comment/uncomment line below to save files
        plt.savefig(f"{path}Plots/{molarity(n)}uM_t{str(trial(n))}_derivs")
        plt.xlabel('Temperature')
    plt.show();

def melting():
    with open(path + 'tM_values.csv', 'w') as file:
        writer = csv.writer(file)
        for n in col:
            names = []
            row = []
            for name in headers[n - 1][:7]:
                s = df_deriv.loc[20:95, name]
                names.append(f"{name}_{molarity(n)}uM_t{trial(n)}")
                row.append(s.idxmin())
                # print(f"{name}_{molarity(n)}uM_t{trial(n)}: {s.idxmin()}")
            writer.writerow(names)
            writer.writerow(row)

# reading the data into panda DataFrames
path = "/Users/carlho/Documents/Shakh Lab/2021_06_21_L82V_hingemuts/"
data = lambda fname : pd.read_excel(path + fname, 
                      sheet_name="FRET", 
                      header=0,
                      index_col=0,
                      usecols="B:CT")
df_amp = data("badkar_2021-06-21 16-31-27_CT010444 - " + 
                " Melt Curve Amplification Results.xlsx")
df_deriv = data("badkar_2021-06-21 16-31-27_CT010444 - " + 
                " Melt Curve Derivative Results.xlsx")

# wanted header names in DataFrames
rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
cols = range(1, 13)
headers = []

for num in cols:
    headers.append([])
    for letter in rows:
        headers[-1].append(letter + str(num))
col = [1, 2, 4, 5]  # columns where the samples are stored
molarity = lambda n : "2" if n < 3 else "4"
trial = lambda n : n if n < 3 else n - 3

# plot();  # uncomment to plot graphs
melting();