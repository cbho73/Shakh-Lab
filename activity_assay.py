from os import lseek, system
import os
from typing import NamedTuple
from numpy.core.fromnumeric import std
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
    count = 0
    for df in dfs:
        plt.figure();
        fig, ax = plt.subplots(2, 4)
        # for row in ax:
        #     for cell in row:
        #         cell.set_xlim([-0.1, 2.2])
        #         cell.set_ylim([0.575, 0.75])
        fig.set_figheight(10)
        fig.set_figwidth(20)

        for col in df.columns:
            c = df.columns.get_loc(col)
            ax[c//4, c % 4].scatter(x=mins, y=df[col]);
            ax[c//4, c % 4].set_title(f"Well {col}: {conc_dict[col]} uM");

        plt.savefig(f"{path}/{order(count)}_{df.index.name}")
        plt.close()
        count = count + 1


def get_v0(change):
    # dictionary that stores the values of point to start at (1-indexed)
    # for key of trial and concentration
    error = {
        # ("P112G T1", 10) : (2, False),
        # ("P112G T1", 25) : (2, False),
        # ("P112G T1", 500) : (3, False),
        # ("V117G T1", 50) : (4, False),
        # ("V117G T1", 200) : (2, False),
        # ("V117G T2", 10) : (2, False),
        # ("L168G T1", 100) : (6, False),
        # ("L168G T1", 500) : (2, False),
        # ("117-164G T1", 100) : (2, False),
        # ("117-164G T1", 400) : (2, False),
        # ("117-164G T1", 500) : (2, False),
        # ("D118N T2", 10) : (3, True)
    }
    fit_params = ["slope", "intercept", "r_value", "p_value", "std_err"]
        # paramters returned from scipy fit to mode of inhibition
    v0_df = pd.DataFrame(columns=conc)  # v0s for the different trials
    regs_data = []  # data from linear regressions for each v0 calculation


    count = 0
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

                # axs.set_xlim([-0.1, 2.2])
                # axs.set_ylim([0.575, 0.75])

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
        fig.savefig(os.path.join(parent_dir, f"v0_test/{order(count)}_{df.index.name}"))
        plt.close()
        regs_data.append(linregs)  # linear regression data appending
        v0_s = pd.Series(v0s, index=conc, name=df.index.name)
        v0_df = v0_df.append(v0_s)  # adding v0s from n trial to df with all v0s
    # outputting data to csv
    v0_df.to_csv(path_or_buf=os.path.join(parent_dir, "v0_values.csv"))
    return v0_df


def fit_v0s(v0_df, dir, order, kMvkI=False, stderr=None, overlay=None):
    params = ["v_max", "k_M", "k_I"]  # parameters for the curves fitted to inhibition
    fit_curves = pd.DataFrame(columns=params)
    path = make_path(parent_dir, f"Activity Plots/{dir}")

    count = 0
    # print(stderr)
    for name, v0_s in v0_df.iterrows():
        # appending to v0s
        # fitting to inhibition model
        try: 
            popt, _pcov = sc.optimize.curve_fit(noncompetitive, 
                                                list(v0_s.index), 
                                                v0_s.to_list(), 
                                                bounds=(0, np.inf), 
                                                verbose=0)
        except RuntimeError:
            print(f"unable to fit for {df.index.name}")

        fitted = pd.Series(data=popt, index=params, name=v0_s.name)
        fit_curves = fit_curves.append(fitted)
    for name, fitted in fit_curves.iterrows():
        # plotting the v0 points and fitted function
        v0_s = v0_df.loc[name]
        plt.figure()
        def fn_range(fn_params):
            pts = []
            for x in np.arange(start=0, stop=conc[-1], step=1):
                pts.append(noncompetitive(x, 
                                          fn_params['v_max'], 
                                          fn_params['k_M'], 
                                          fn_params['k_I']))
            return pts
        plt.plot(fn_range(fitted), label=name)
        plt.scatter(x=conc, y=v0_s)
        if stderr is not None:
            plt.errorbar(x=conc, y=v0_s, fmt='o', yerr=stderr.loc[name])
        plt.savefig(f"{path}/{order(count)}_{name}")
        if overlay is not None and fitted.loc['k_I'] > fit_curves.loc[overlay, 'k_I']:
            plt.plot(fn_range(fit_curves.loc[overlay]), label=overlay)
            plt.legend()
            overlay_path = make_path(path, 'overlay')
            plt.savefig(f"{overlay_path}/{order(count)}_{name}")
        plt.close()
        count = count + 1
    
    # plotting k_M vs k_I
    if kMvkI:
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel(r"$K_M$")
        ax.set_ylabel(r"$K_I$")
        for name, row in fit_curves.iterrows():
            ax.scatter(x=row['k_M'], y=row['k_I'], label=name)
        ax.legend()
        fig.savefig(f"{parent_dir}/kM_vs_kI")
        plt.close()
    fit_curves.to_csv(path_or_buf=os.path.join(path, "curve_params.csv"))
    return fit_curves

def average_v0(v0_df):
    avg_df = pd.DataFrame()  # v0s for the different trials
    stdev_df = pd.DataFrame()  # standard deviation for the trials
    for t in range(9):
        name = v0_df.index[t*ntrials].strip().split()[0]
        if t != 8:
            avg = v0_df[t*ntrials:t*ntrials+ntrials].mean(axis=0)
            stdev = v0_df[t*ntrials:t*ntrials+ntrials].std(axis=0)
            stdev.name = name
            avg.name = name
        else:
            avg = v0_df[t*ntrials:t*ntrials+2].mean(axis=0)
            stdev = v0_df[t*ntrials:t*ntrials+2].std(axis=0)
            stdev.name = name
            avg.name = name
        avg_df = avg_df.append(avg)
        stdev_df = stdev_df.append(stdev)
    order = lambda count : ordering[count]
    fit_v0s(avg_df, 
            "avg v0 plots", 
            order, 
            kMvkI=True, 
            stderr=stdev_df, 
            overlay='L82V')

# function that describes mode of enzymatics to fit v_0 data
def noncompetitive(amp, v_max, k_M, k_I):
    return (v_max * amp) / (k_M + amp * (1 + amp/k_I))

# PARAMETERS TO CHANGE
lab_dir = "/Users/carlho/Documents/Shakh Lab/"
fname = "L82V_muts_activity_2.xlsx"
run = 2  # what run of activity assay is this?
num_tables = 35  # number of sheets in excel
table_len = 25  # length of table data entries
spacers = 2  # length of non-data parts of table, header, spacing, etc
headers = [(table_len + spacers) * n + 2 for n in range(num_tables)]
ntrials = 4
ordering = [5, 6, 9, 1, 2, 3, 4, 7, 8]  # trial numbering (because I went out of order)
order = lambda count : ordering[int(count / ntrials)]
trials = [2 if n > 7 else 4 for n in range(9)]
trial_order = dict(zip(ordering, trials))

label = ["A", "B", "C", "D", "E", "F", "G", "H"] # labels for wells
conc = [2.5, 5, 10, 25, 50, 100, 250, 500]     # concentrations
conc_dict = dict(zip(label, conc))               

dfs = []
parent_dir = make_path(lab_dir, f"activity_run{run}")
make_path(parent_dir, "Activity Plots")
for header in headers:
    df = pd.read_excel(os.path.join(lab_dir, fname), 
                             sheet_name="Sheet1", 
                             header=header,
                             index_col=0, 
                             usecols="A:I",
                             nrows=table_len)
    df.columns = label
    dfs.append(df)

mins = pd.Series(dfs[0].index).apply(to_minutes) # convert time to mins

redo = True
v0_path = os.path.join(parent_dir, "v0_values.csv")
if os.path.exists(v0_path) and not redo:
    v0_df = pd.read_csv(v0_path, 
                             header=0,
                             index_col=0, 
                             usecols=[n for n in range(9)])
    # print(v0_df)
else:
    v0_df = get_v0(8)

# plot_traces(parent_dir, "Activity Plots/Traces")
# fit_v0s(v0_df, "v0 plots", order)
average_v0(v0_df)