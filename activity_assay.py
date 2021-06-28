from os import lseek
import os
from numpy.lib.twodim_base import tri
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import scipy as sc
from scipy import stats
import numpy as np

def make_path(parent, child):
    path = os.path.join(parent, child)
    if not os.path.isdir(path):
        os.mkdir(path)
    return path

parent_dir = "/Users/carlho/Documents/Shakh Lab/"
data = lambda fname, sname : pd.read_excel(parent_dir + fname, 
                                           sheet_name=sname, 
                                           header=0,
                                           index_col=0,
                                           usecols="A:I")
dfs = []
for n in range(16):
    dfs.append(data("L82V_muts_activity_py.xlsx", f"Sheet{n + 1}"))

def to_minutes(t):
    o = dt.time(0, 0, 0)
    date = dt.date(1,1,1)
    combine = lambda time : dt.datetime.combine(date, time)
    t_delta = combine(t) - combine(o)
    return t_delta.seconds/60

mins = pd.Series(dfs[0].index).apply(to_minutes) # convert time to mins
label = ["A", "B", "C", "D", "E", "F", "G", "H"] # labels for wells
conc = [10, 25, 50, 100, 200, 300, 400, 500] # concentrations
conc_dict = dict(zip(label, conc))

def plot_traces():
    path = make_path(parent_dir, "Activity Plots/Traces")
    count = 2
    for df in dfs:
        plt.figure();
        # axes = plt.gca()
        # axes.set_xlim([xmin,xmax])
        fig, ax = plt.subplots(2, 4)
        for row in ax:
            for cell in row:
                cell.set_ylim([0.575, 0.75])
        fig.set_figheight(10)
        fig.set_figwidth(20)
        c = 0
        for col in df.columns:
            ax[c//4, c % 4].scatter(x=mins, y=df[col]);
            ax[c//4, c % 4].set_title(f"Well {col}: {conc_dict[col]} uM");
            c = c + 1
        plt.savefig(f"{path}/{int(count / 2)}_{df.index.name}")
        count = count + 1


def get_v0(change):
    # dictionary that stores the values of point to start at (1-indexed)
    # for key of trial and concentration
    error = {
        ("P112G T1", 10) : (2, False),
        ("P112G T1", 25) : (2, False),
        ("P112G T1", 500) : (3, False),
        ("V117G T1", 50) : (4, False),
        ("V117G T1", 200) : (2, False),
        ("V117G T2", 10) : (2, False),
        ("L168G T1", 100) : (6, False),
        ("L168G T1", 500) : (2, False),
        ("117-164G T1", 100) : (2, False),
        ("117-164G T1", 400) : (2, False),
        ("117-164G T1", 500) : (2, False),
        ("D118N T2", 10) : (3, True)
    }
    params = ["v_max", "k_M", "k_I"]  # parameters for the curves fitted to inhibition
    fit_params = ["slope", "intercept", "r_value", "p_value", "std_err"]
        # paramters returned from scipy fit to mode of inhibition

    v0_df = pd.DataFrame(columns=conc)  # v0s for the different trials
    regs_data = []  # data from linear regressions for each v0 calculation
    fit_curves = pd.DataFrame(columns=params)

    count = 2
    for df in dfs:
        linregs = pd.DataFrame(columns=fit_params)
        plt.figure()
        v0s = []
        for col in df.columns:
            trace = pd.Series(df[col]).to_list()
            def line_slope(start, delta):
                t = mins[start : start + delta]
                c = trace[start : start + delta]
                lin_eq = list(stats.linregress(t, c))
                # print(lin_eq)
                df_line = pd.DataFrame([lin_eq], index=[df.index.name], columns=fit_params)
                linregs.append(df_line)
                return lin_eq[0]
            try:
                p, b = error[df.index.name, col]
                if not b:
                    v0s.append(-line_slope(p - 1, p - 1 + change))
                else:
                    slope = line_slope(p - 1, 5)
                    diff = trace[p - 2] - (slope * (mins[p - 2] - mins[p - 1]) + mins[p - 1])
                    for i in range(p - 1, len(trace)):
                        trace[i] = trace[i] + diff
                    v0s.append(-line_slope(0, change))
            except KeyError:
                v0s.append(-line_slope(0, change))
            
        regs_data.append(linregs)  # linear regression data appending

        # appending to v0s
        v0_s = pd.Series(v0s, index=conc, name=df.index.name)
        v0_df = v0_df.append(v0_s)  # adding v0s from n trial to df with all v0s

        # fitting to inhibition model
        popt, _pcov = sc.optimize.curve_fit(noncompetitive, 
                                    list(v0_s.index), 
                                    v0_s.to_list())
        fitted = pd.Series(data=popt, index=params, name=v0_s.name)
        fit_curves = fit_curves.append(fitted)

        # plotting the points and function
        plt.figure()
        pts = []
        for x in np.arange(start=0, stop=conc[-1], step=1):
            pts.append(noncompetitive(x, fitted['v_max'], fitted['k_M'], fitted['k_I']))
        plt.plot(pts)
        plt.scatter(x=conc, y=v0s)
        path = make_path(parent_dir, "Activity Plots/v0 plots")
        plt.savefig(f"{path}/{int(count/2)}_{df.index.name}")
        count = count + 1
    v0_df.to_csv(path_or_buf=os.path.join(parent_dir, "v0_values.csv"))
    fit_curves.to_csv(path_or_buf=os.path.join(parent_dir, "curve_params.csv"))

# function that describes mode of enzymatics to fit v_0 data
def noncompetitive(amp, v_max, k_M, k_I):
    return (v_max * amp) / (k_M + amp * (1 + amp/k_I))
        
path = make_path(parent_dir, "Activity Plots/")

plot_traces()
# v0_df = get_v0(8)