from os import lseek, system
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

def to_minutes(t):
    o = dt.time(0, 0, 0)
    date = dt.date(1,1,1)
    combine = lambda time : dt.datetime.combine(date, time)
    t_delta = combine(t) - combine(o)
    return t_delta.seconds/60

def plot_traces(parent_dir, new_dir):
    path = make_path(parent_dir, new_dir)
    count = 2
    for df in dfs:
        plt.figure();
        # axes = plt.gca()
        # axes.set_xlim([xmin,xmax])
        fig, ax = plt.subplots(2, 4)
        for row in ax:
            for cell in row:
                cell.set_xlim([-0.1, 2.2])
                cell.set_ylim([0.575, 0.75])
        fig.set_figheight(10)
        fig.set_figwidth(20)

        for col in df.columns:
            c = df.columns.get_loc(col)
            ax[c//4, c % 4].scatter(x=mins, y=df[col]);
            ax[c//4, c % 4].set_title(f"Well {col}: {conc_dict[col]} uM");

        plt.savefig(f"{path}/{int(count / 2)}_{df.index.name}")
        plt.close()
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
        v0s = []

        plt.figure()
        fig, ax = plt.subplots(2, 4)
        fig.set_figheight(10)
        fig.set_figwidth(20)

        for col in df.columns:
            trace = pd.Series(df[col]).to_list()
            def line_slope(start, delta):
                t = mins[start : start + delta]
                y_conc = trace[start : start + delta]
                lin_eq = list(stats.linregress(t, y_conc))
                df_line = pd.DataFrame([lin_eq], index=[df.index.name], columns=fit_params)
                linregs.append(df_line)

                # testing (outputting points and line)
                c = df.columns.get_loc(col)
                axs = ax[c//4, c % 4]
                axs.scatter(x=t, y=y_conc);
                axs.set_title(f"Well {col}: {conc_dict[col]} uM");
                axs.set_xlim([-0.1, 2.2])
                axs.set_ylim([0.575, 0.75])

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
        
        make_path(parent_dir, "v0_test")
        fig.savefig(os.path.join(parent_dir, f"v0_test/{int(count / 2)}_{df.index.name}"))
        plt.close()

        regs_data.append(linregs)  # linear regression data appending

        # appending to v0s
        v0_s = pd.Series(v0s, index=conc, name=df.index.name)
        v0_df = v0_df.append(v0_s)  # adding v0s from n trial to df with all v0s

        # fitting to inhibition model
        popt, _pcov = sc.optimize.curve_fit(noncompetitive, 
                                            list(v0_s.index), 
                                            v0_s.to_list(), 
                                            bounds=(0, np.inf), 
                                            verbose=0)

        fitted = pd.Series(data=popt, index=params, name=v0_s.name)
        fit_curves = fit_curves.append(fitted)

        # plotting the v0 points and fitted function
        plt.figure()
        pts = []
        for x in np.arange(start=0, stop=conc[-1], step=1):
            pts.append(noncompetitive(x, fitted['v_max'], fitted['k_M'], fitted['k_I']))
        plt.plot(pts)
        plt.scatter(x=conc, y=v0s)
        path = make_path(parent_dir, "Activity Plots/v0 plots")
        plt.savefig(f"{path}/{int(count/2)}_{df.index.name}")
        plt.close()
        count = count + 1
    
    # plotting k_M vs k_I
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel(r"$K_M$")
    ax.set_ylabel(r"$K_I$")
    for name, row in fit_curves.iterrows():
        ax.scatter(x=row['k_M'], y=row['k_I'], label=name)
    ax.legend()
    fig.savefig(f"{parent_dir}/kM_vs_kI")
    plt.close()

    # outputting data to csv
    v0_df.to_csv(path_or_buf=os.path.join(parent_dir, "v0_values.csv"))
    fit_curves.to_csv(path_or_buf=os.path.join(parent_dir, "curve_params.csv"))

# function that describes mode of enzymatics to fit v_0 data
def noncompetitive(amp, v_max, k_M, k_I):
    return (v_max * amp) / (k_M + amp * (1 + amp/k_I))

parent_dir = "/Users/carlho/Documents/Shakh Lab/"
make_path(parent_dir, "Activity Plots")
fname = "L82V_muts_activity_py.xlsx"
data = lambda sname : pd.read_excel(parent_dir + fname, 
                                           sheet_name=sname, 
                                           header=0,
                                           index_col=0,
                                           usecols="A:I")

num_sheets = 16  # number of sheets in excel
dfs = []
for n in range(num_sheets):
    dfs.append(data(f"Sheet{n + 1}"))

mins = pd.Series(dfs[0].index).apply(to_minutes) # convert time to mins
label = ["A", "B", "C", "D", "E", "F", "G", "H"] # labels for wells
conc = [10, 25, 50, 100, 200, 300, 400, 500]     # concentrations
conc_dict = dict(zip(label, conc))               

# plot_traces(parent_dir, "Activity Plots/Traces")
v0_df = get_v0(8)