import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import pandas as pd
from opt import mystyle


# Make plots nicer!
plt.style.use("seaborn-colorblind")
plt.rcParams.update(mystyle())

cells = ["cck", "olm", "vipcck", "vipcr"]

# Variable to store the results
cols = ["activity", "spikes", "cell"]
df = pd.DataFrame(columns=cols, index=range(len(cells)*10))
counter = 0
replications = 10
dirname = "data_interneurons_synapses/"

df2 = pd.DataFrame(columns=["peak", "cell"], index=range(len(cells)*10*200))

counter2 = 0
for cell in cells:
    
    spikes = []
    silent = []
    active = []
    with open(dirname + cell + '.pkl', 'rb') as handle:
        data = pickle.load(handle)
    
    for j in range(10):
        data_ = data[f'Replication_{j}_voltage_soma']
        L = len(data_)
        maxs = [np.max(vol) for vol in data_ if np.max(vol)>30]
        sils = [np.max(vol) for vol in data_ if np.max(vol)<1e-1]
        
        silent.append(len(sils))
        active.append((200-len(sils)))
        spikes.append(len(maxs)/active[-1])

        df.loc[counter].activity = active[-1]/200*100
        df.loc[counter].spikes = spikes[-1]*100
        df.loc[counter].cell = cell
        counter += 1
        
        maxs2 = [np.max(vol) for vol in data_]
        for jj in range(200):
            df2.loc[counter2].peak = maxs2[jj]
            df2.loc[counter2].cell = cell
            counter2 += 1
        

# df_graph = df2["peak"][df2["cell"]=="vipcck"]
# test = np.array(df_graph)
# test2 = test.reshape(10, 200).T


plt.figure()
plt.subplot(121)
sns.barplot(data=df, y="activity", x="cell", ci="sd", capsize=.2,
            palette="Set2")
plt.ylim([0, 100])
plt.ylabel("Percent (%)")
plt.xlabel("cell type")
plt.title("LEC induced activity")
plt.subplot(122)
sns.barplot(data=df, y="spikes", x="cell", ci="sd", capsize=.2,
            palette="Set2")
plt.ylabel("Percent (%)")
plt.xlabel("cell type")
plt.title("LEC induced spikes")
plt.tight_layout()
sns.color_palette("Set2")
plt.show()