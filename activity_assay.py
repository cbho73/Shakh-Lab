from os import lseek
import os
from numpy.lib.twodim_base import tri
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import scipy

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

time_mins = pd.Series(dfs[0].index).apply(to_minutes)
conc_dict = {
    "A" : 10,
    "B" : 25,
    "C" : 50,
    "D" : 100,
    "E" : 200,
    "F" : 300,
    "G" : 400,
    "H" : 500
}
conc = [10, 25, 50, 100, 200, 300, 400, 500]

def plot_traces():
    trace_dir = "Activity Plots/Traces"
    path = os.path.join(parent_dir, trace_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    count = 2
    for df in dfs:
        plt.figure();
        fig, ax = plt.subplots(2, 4)
        fig.set_figheight(10)
        fig.set_figwidth(20)
        c = 0
        for col in df.columns:
            ax[c//4, c % 4].scatter(x=time_mins, y=df[col]);
            ax[c//4, c % 4].set_title(conc_dict[col]);
            c = c + 1
        plt.savefig(f"{path}/{int(count / 2)}_{df.index.name}")
        count = count + 1

def get_v0():
    for df in dfs:
        plt.figure()
        for col in df.columns:
            trace = pd.Series(df[col]).to_list()
            

dir = "Activity Plots"
path = os.path.join(parent_dir, dir)
if not os.path.isdir(path):
    os.mkdir(path)


plot_traces()

