import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from activity_assay import make_path

path = "/Users/carlho/Documents/Shakh Lab/FLUORIMETRY"
names = [f"L82V_muts_{n}" for n in range(1, 7)] + ["L82V", "CONTROL"]
dilutions = [6, 5, 4, 3, 2, 1, 0, 0]
names_dict = dict(zip(names, dilutions))
dfs = []
for name in names:
    if 'L82V' in name:
        if 'muts' in name:
            col = 8
        else:
            col = 6
    else:
        col = 2
    df = pd.read_csv(os.path.join(path, name + '.csv'), 
                     header=[0, 1], 
                     nrows=146)

    # df.columns = pd.MultiIndex.from_tuples([tuple('' if y.startswith('Unnamed:') else y for y in x) for x in map(list, df.columns.tolist())])
    df.index.name = name
    # print(df)
    dfs.append(df)
conc = lambda n, init : int((init * 0.5 * (0.4 ** n)) * 1000) / 1000
init_conc = {
    '117-164G' : 198,
    '164G' : 247,
    '117G' : 215,
    'L82V' : 62.7,
    'WT' : 1064.6,
}
prot_names = ['117-164G', '164G', '117G', 'WT']
colors = cm.rainbow(np.linspace(0, 1, len([n for n in range(6)])))
data_dir = make_path(path, "graphs")
for name in prot_names:
    plt.figure()
    plt.title(name)
    for i, c in zip([n for n in range(6)], colors):
        # print("trial: " + str(i))
        # print(dfs[i])
        plt.scatter((dfs[i])[name, 'Wavelength (nm)'], 
                    (dfs[i])[name, 'Intensity (a.u.)'], 
                    color=c,
                    label=conc(dilutions[i], init_conc[name]))
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='xx-small')
    plt.savefig(fname=os.path.join(data_dir, name),
                bbox_extra_artists=(lgd,), 
                bbox_inches='tight')
    plt.close()
colors = cm.rainbow(np.linspace(0, 1, 4))
plt.figure()
plt.title('L82V')
for i, c in zip(np.arange(1, 5), colors):
    plt.scatter((dfs[6])[str(i), 'Wavelength (nm)'], 
                (dfs[6])[str(i), 'Intensity (a.u.)'], 
                color=c,
                label=conc(dilutions[i], init_conc['L82V']))
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='xx-small')
    plt.savefig(fname=os.path.join(data_dir, 'L82V'),
                bbox_extra_artists=(lgd,), 
                bbox_inches='tight')
plt.close()