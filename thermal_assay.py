import pandas as pd
import matplotlib.pyplot as plt

rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
cols = range(1, 13)
headers = []

for num in cols:
    headers.append([])
    for letter in rows:
        headers[-1].append(letter + str(num))
print(len(headers))

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

col = [1, 2, 4, 5]  # columns where the samples are stored
for n in col:
    df1 = df_amp.loc[20:95, headers[n - 1][:7]]
    df2 = df_deriv.loc[20:95, headers[n - 1][:7]]
    molarity = "2" if n < 3 else "4"
    trial = n if n < 3 else n - 3
    
    df1.plot();
    plt.savefig(f"{path}Plots/{molarity}uM_t{str(trial)}")
    df2.plot();
    plt.savefig(f"{path}Plots/{molarity}uM_t{str(trial)}_derivs")
    plt.xlabel('Temperature')
plt.show();