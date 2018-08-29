import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from natsort import natsorted

def pretty_plotter(df, proteins, labels, concs, markers, xlabel, title, output, spacing=4, savefig=False):

    summary_df = pd.DataFrame(columns=df.columns)
    summary_df["biological reps"] = None

    for i in range(len(proteins)):
        for j in range(len(concs)):

            x = df.loc[(df["protein"] == proteins[i]) & (df["conc"] == concs[j])]

            if concs[j] in x.conc.values:
                summary_df.at[len(concs)*i+j, "average"] = np.mean(x.average)
                summary_df.at[len(concs)*i+j, "stderr"] = np.std(x.average)/np.sqrt(len(x.average))
                summary_df.at[len(concs)*i+j, "biological reps"] = int(len(x.average))
                summary_df.at[len(concs)*i+j, "color"] = x.color.values[0]
                summary_df.at[len(concs)*i+j, "conc"] = x.conc.values[0]
                summary_df.at[len(concs)*i+j, "label"] = x.label.values[0]
                summary_df.at[len(concs)*i+j, "protein"] = x.protein.values[0]
            else:
                continue

    summary_df.reset_index(inplace=True, drop=True)

    for i in range(len(labels)):
        plt.errorbar(summary_df.conc[spacing*i:spacing*i+spacing], summary_df.average[spacing*i:spacing*i+spacing],
                     yerr=summary_df.stderr[spacing*i:spacing*i+spacing], marker=markers[i], mfc=summary_df.color[spacing*i+1],
                     mec="black", ms=8, ecolor="black", capsize=3, fmt="-k")

    plt.legend(labels, fontsize=10, loc="center right", markerscale=0.75, bbox_to_anchor=(1.2,0.5),fancybox=True, shadow=True)
    plt.plot([2,12], [100, 100], "k--")
    plt.xticks(fontsize=14)
    plt.xlabel(xlabel, fontsize=16)
    plt.yticks(fontsize=14)
    plt.ylabel("% of untreated S. epi \n growth at 12 hours", fontsize=16)
    plt.ylim(0,120)
    plt.title(title, fontsize=18)
    plt.tight_layout()
    if savefig == True:
        plt.savefig(output, format='svg', dpi=300)
    None

    return summary_df
