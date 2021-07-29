from activity_assay import make_path
from os import lseek
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import importlib
import numpy as np
import scipy as sc
import scipy.signal as sg
import os
importlib.import_module("activity_assay")

# plot and save graphs
def plot():
    # print(df_amp.head())
    make_path(path, "Plots")
    make_path(path, "Plots/derivs")
    make_path(path, "Plots/amp")
    colors = cm.rainbow(np.linspace(0, 1, len(col)))
    for n in col:
        df1 = df_amp.loc[20:95, headers[n - 1]]
        df2 = df_deriv.loc[20:95, headers[n - 1]]
        
        plt.figure();
        df1.plot(cmap='magma');
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='xx-small')
        plt.savefig(f"{path}Plots/amp/{n}_{prot_substrate(n)}_t{str(trial(n))}.png", 
                    bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
        plt.xlabel('Temperature')
        plt.close();

        plt.figure();
        df2.plot(cmap='rainbow');
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='xx-small')
        plt.savefig(f"{path}Plots/derivs/{n}_{prot_substrate(n)}_t{str(trial(n))}_derivs.png", 
                    bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
        plt.xlabel('Temperature')
        plt.close()
    # plt.show();

def melting():
    AMP_ATP = lambda n : AMP_conc if n < ATP_start else ATP_conc
    AMP_tM = []
    ATP_tM = []
    AMP_names = []
    ATP_names = []
    newpath = make_path(path, "tM_Plots")
    params = ['top', 'bottom', 'hillslope', 'x_0']
    lst_params = []
    for n in col:
        name = f"{n}_{prot_substrate(n)}_t{trial(n)}"
        row = []
        ignored = ignore.get(n, [])
        count = 1
        for i in headers[n - 1]:
            def check(strng):
                for i in [7, 8, 17, 18]:
                    if str(i) in strng:
                        return True
                return False
            
            if check(i):
                print("success")
                s = df_deriv.loc[54:65, i]
                print(s)
            else:
                s = df_deriv.loc[45:65, i]
            try:
                ignored.index(count)
                print(f"ignored point {i} in column {n}")
            except ValueError:
                # lst = list(s.map(lambda n : n * -1))
                # print(lst)
                # peaks, _data = sg.find_peaks(lst)
                # print(list(peaks))
                # print(s)
                # mins = s.nsmallest(2).sort_index()
                # if abs(s.index[0] - s.index[1]) > 0.5:
                #     print(mins)
                #     print(abs(s.index[0] - s.index[1]))
                #     row.append(mins.iloc[0])
                # else:
                #     row.append(s.idxmin())
                row.append(s.idxmin())
            count = count + 1
            # print(f"{name}_{molarity(n)}uM_t{trial(n)}: {s.idxmin()}")

        conc = AMP_ATP(n).copy()
        try: 
            for i in ignored:
                conc.remove(AMP_ATP(n)[i - 1])
            # plt.figure()
            # plt.title(name)
            # plt.xscale("log")
            # plt.scatter(conc, row)
            # plt.show()
            # plt.close()
            if len(ignored) == 0:
                popt, _pcov = sc.optimize.curve_fit(sigmoid, 
                                                    conc, 
                                                    row,
                                                    bounds=(0, np.inf), 
                                                    verbose=0)
            else:
                popt, _pcov = sc.optimize.curve_fit(sigmoid, 
                                                    conc, 
                                                    row)
            lst_params.append(popt)
        except RuntimeError:
            popt = [None] * 4
            lst_params.append(popt)
            print(f"unable to fit for {name}")

        def fn_range(fn_params):
            pts = []
            iter = np.arange(start=-5, stop=2, step=0.25, dtype=np.float128)
            molarity = []
            for x in iter:
                # print(10**x)
                molarity.append(10**x)
                pts.append(sigmoid(10**x, 
                                    fn_params['top'], 
                                    fn_params['bottom'], 
                                    fn_params['hillslope'],
                                    fn_params['x_0']))
            # exp = lambda x : 10**x
            return pd.Series(data=pts, index=molarity)
        
        temp = row.copy()
        for i in ignored:
            temp.insert(i - 1, None)
        if n < ATP_start :
            AMP_tM.append(temp)
            AMP_names.append(name)
        else :
            ATP_tM.append(temp)
            ATP_names.append(name)

        plt.figure()
        plt.xscale("log")
        if popt[0] is not None:
            fit = pd.Series(data=popt, index=params, name=name)
            fn = fn_range(fit)
            plt.plot(fn.index, list(fn))
        plt.scatter(conc, row)
        # print(f"{newpath}/{n}_{prot_substrate(n)}_t{trial(n)}")
        plt.savefig(f"{newpath}/{n}_{prot_substrate(n)}_t{trial(n)}")
        plt.close()
    AMP_df = pd.DataFrame(AMP_tM, index=AMP_names, columns=AMP_conc)
    ATP_df = pd.DataFrame(ATP_tM, index=ATP_names, columns=ATP_conc)
    AMP_df.to_csv(path_or_buf=os.path.join(newpath, "AMP_tM.csv"))
    ATP_df.to_csv(path_or_buf=os.path.join(newpath, "ATP_tM.csv"))
    # print(AMP_df)
    # print(ATP_df)
    fit_params = pd.DataFrame(lst_params, index=(AMP_names + ATP_names), columns=params)
    fit_params.to_csv(path_or_buf=os.path.join(newpath, "curve_params.csv"))

def sigmoid(conc, top, bottom, hillslope, x_0):
    return bottom + (top-bottom)/(1+(x_0/conc)**hillslope)

def collate(lst):
    temp = []
    half = int(len(lst) / 2)
    if len(lst) % 2 != 0:
        raise ValueError('cannot use odd-numbered lists for collate.')
    for i in range(half):
        temp.append(lst[i * 2])
    for i in range(half):
        temp.append(lst[(i * 2) + 1])
    return temp

def cat_str(lst, str):
    temp = []
    for elt in lst:
        temp.append(elt + str)
    return temp

def alternate(df, cols):
    all_headers = []
    ord_headers = []
    for i in cols:
        r = cat_str(rows, str(i))
        ord_headers = ord_headers + r
        all_headers = all_headers + collate(r)
    # print(all_headers)
    # print(df)
    df = df[all_headers]
    # print(df)
    df.columns = ord_headers
    # print(df)
    return df
# reading the data into DataFrames
path = "/Users/carlho/Documents/Shakh Lab/20210728_TSA/"
data = lambda fname : pd.read_excel(path + fname, 
                      sheet_name="FRET", 
                      header=0,
                      index_col=0,
                      usecols="B:NV")
df_amp = data("badkar_2021-07-28 04-27-42_CT010451 -  Melt Curve Amplification Results.xlsx")
df_deriv = data("badkar_2021-07-28 04-27-42_CT010451 -  Melt Curve Derivative Results.xlsx")

ntrials = 20 # number of trials
# wanted header names in DataFrames
rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
cols = range(1, 25)
headers = []
ATP_start = (ntrials / 2) + 1 # column number where ATP starts
for num in cols:
    headers.append([])
    for letter in rows:
        headers[-1].append(letter + str(num))
prot_key = {
    1 : "L82V-117G",
    2 : "L82V-164G",
    3 : "L82V-117-164G", 
    4 : "L82V", 
    5 : "WT", 
    0 : "WT"
}
trial_spacing = {
    'close' : (lambda n : ((n + 1) % 2) + 1, 
               lambda n : prot_key[int(((n % (ntrials/2)) + 1) / 2)] + ("_AMP" if n < ATP_start else "_ATP"),
               [n for n in range(1, ntrials + 1)]),
    'e_o' : (lambda n : (int((n - 1) / 2) % 2) + 1, 
             lambda n : prot_key[int(((((n + 1) / 2) % (ntrials/2)) + 1) / 2)] + ("_AMP" if n < ATP_start else "_ATP"),
             [2 * n + 1 for n in range(ntrials)])
}
trial, prot_substrate, col = trial_spacing['close']

# # FOR THE NIGHTMARE THAT IS ALTERNATE CONCENTRATIONS
df_amp = alternate(df_amp, cols)
df_deriv = alternate(df_deriv, cols)

# well concentrations for AMP and ATP
AMP_init = 72.5
ATP_init = 119
AMP_conc = []
ATP_conc = []
for n in range(16):
    AMP_conc.append(AMP_init * 0.5 * (0.4 ** n))
    ATP_conc.append(ATP_init * 0.5 * (0.4 ** n))

# molarity = lambda n : "2" if n < 3 else "4"
# trial = lambda n : n if n < 3 else n - 3

# points to ignore
ignore = {
    1 : [1, 2],
    5 : [1],
    7 : [1, 2],
    9 : [1, 2],
    11 : [2],
    13 : [2],
    15 : [2],
    17 : [1, 2, 3],
    21 : [1, 2],
}
ignore = {
    2 : [2],
    5 : [1, 16],
    6 : [1, 2],
    8 : [4],
    9 : [1, 15, 16],
    10 : [1, 15, 16],
    11 : [16],
    12 : [16],
    13 : [1, 2, 3, 13, 15, 16],
    14 : [1, 2, 3, 4, 15, 16],
    15 : [1, 16],
    16 : [1, 16],
}
ignore = {}
# print(df_amp.index)

plot();  
# uncomment to plot graphs
melting();
