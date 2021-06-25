from os import lseek
import os
from numpy.lib.twodim_base import tri
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
from scipy import stats

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
                                    usecols="A:I"
                                    )
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
        fig, ax = plt.subplots(2, 4)
        fig.set_figheight(10)
        fig.set_figwidth(20)
        c = 0
        for col in df.columns:
            ax[c//4, c % 4].scatter(x=mins, y=df[col]);
            ax[c//4, c % 4].set_title(f"Well {col}: {conc_dict[col]} uM");
            c = c + 1
        plt.savefig(f"{path}/{int(count / 2)}_{df.index.name}")
        count = count + 1

def get_v0():
    # dictionary that stores the values of point to start at (1-indexed)
    # for key of trial and concentration
    start = {
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
    count = 2
    for df in dfs:
        plt.figure()
        v0s = []
        for col in df.columns:
            trace = pd.Series(df[col]).to_list()
            def line_slope(start, delta):
                t = mins[start : start + delta]
                c = trace[start : start + delta]
                slope, _intercept, _r_value, _p_value, _std_err = stats.linregress(t, c)
                return slope
            try:
                p, b = start[df.index.name, col]
                if not b:
                    v0s.append(-line_slope(p - 1, 9))
                else:
                    slope = line_slope(p - 1, 5)
                    diff = trace[p - 2] - (slope * (mins[p - 2] - mins[p - 1]) + mins[p - 1])
                    for i in range(p - 1, len(trace)):
                        trace[i] = trace[i] + diff
                    plt.figure()
                    plt.scatter(x=mins, y=trace)
                    plt.show()
                    v0s.append(-line_slope(0, 9))
            except KeyError:
                v0s.append(-line_slope(0, 9))
        plt.figure()
        plt.scatter(x=conc, y=v0s)
        path = make_path(parent_dir, "Activity Plots/v0 plots")
        plt.savefig(f"{path}/{int(count/2)}_{df.index.name}")
        count = count + 1
        
path = make_path(parent_dir, "Activity Plots/")

plot_traces()
get_v0()