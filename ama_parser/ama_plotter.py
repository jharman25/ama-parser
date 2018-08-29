import pandas as pd
import matplotlib.pyplot as plt
from natsort import natsorted

def ama_aggregate_sorter(dataframes):

    df = pd.concat(dataframes)
    df.conc = df.conc.astype('category')
    df.conc.cat.reorder_categories(natsorted(set(df.conc)), inplace=True, ordered=True)
    df = df.sort_values(by=["color", "protein", "conc"])
    df = df.reset_index()
    df = df.drop(df.columns[0], axis=1)

    return df


def ama_plotter(df, plot_title, colors, figsize=(7,7)):

    fig = plt.figure(figsize=figsize)
    plt.bar(range(len(df)), df.average, align='center', yerr=df.stderr, color=colors, linewidth=1, edgecolor="black")
    plt.xticks(range(len(df)), df.label, fontsize=20)
    plt.yticks(fontsize=20)

    plt.plot([0,int(len(df))], [100, 100], "k--")

    plt.ylabel("% growth at 12 hours", fontsize=20)
    plt.title(plot_title, fontsize=25)
    plt.ylim(0, 125)
    plt.tight_layout()
    plt.show()
