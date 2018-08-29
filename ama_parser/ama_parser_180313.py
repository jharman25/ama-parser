import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from natsort import natsorted

def ama_parser_180313(file_name, discard, blanks, proteins, concentrations, time_point, plot_title,
               colors, skiprows=2, skipfooter=14, time_int = 0.25, figsize=(7,7), output = "output_file.png",
               empties = True, averaging = False, growth_curves = False,
               plots = True, save = False):

    #read in data with proper encoding, skipping metadata in first two and last four rows
    df = pd.read_table(file_name, encoding = "utf-16", skiprows = skiprows, skipfooter = skipfooter, engine = 'python')

    #discard wells we don't care about, as specified in list of strings "discard"
    #also automatically discard pre-specified empty outer wells unless empties == False

    df = df.dropna(axis = 1)

    if empties == True:

        junk_and_water_wells = ["Temperature(¡C)", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",
               "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "B1", "C1", "D1", "E1",
               "F1", "G1", "B12", "C12", "D12", "E12", "F12", "G12"]

        df = df.drop(junk_and_water_wells, axis = 1)
        df = df.drop(discard, axis = 1)

    #Rename blank columns - blanks specified by user
    df = df.rename(columns = {blanks[0]:"Blank 1", blanks[1]: "Blank 2", blanks[2]:"Blank 3"})

    #Take average of 3 blanks and subtract from all data
    blank = (df["Blank 1"] + df["Blank 2"] + df["Blank 3"])/3

    d2 = []

    for i in range(len(df.columns)):

        if i+1 < len(df.columns):
            d = df.iloc[:,i+1] - blank
            d2.append(d)

    df = pd.DataFrame(d2).T

    #rename time column in dataframe
    time = []

    for i in range(len(df)):
        if i == 0:
            time.append(i)
        else:
            time.append(time[i-1] + time_int)

    df["Time (hours)"] = time

    data_labels = []

    for i in range(len(proteins)):
        label = concentrations[i] + " uM " + proteins[i]
        data_labels.append(label)

    #optional plotting of triplicate-averaged growth curves when growth_curves == True
    if growth_curves == True:

        for i in range(len(data_labels)):

            means = []
            stderrs = []

            for j in range(len(df)):
                x = np.mean(df[df.columns[i*3+3:i*3+6]].values[j])
                y = np.std(df[df.columns[i*3+3:i*3+6]].values[j])/len(df[df.columns[i*3+3:i*3+6]].values[j])
                means.append(x)
                stderrs.append(y)

            plt.errorbar(df["Time (hours)"], means, yerr=stderrs, fmt=".")
            plt.errorbar(df["Time (hours)"], means, yerr=stderrs, fmt=".")
            if i == 0:
                plt.title("Wildtype growth")
            else:
                plt.title(data_labels[i-1])
            plt.title(data_labels[i])
            plt.xlabel("Time (hours)")
            plt.ylabel("OD600")
            plt.ylim(0,1.5)
            plt.show()

    #default boxplots of all data (leaves out blank)

    if plots == True:

        aves = []

        t = df.index[df["Time (hours)"] == time_point].values

        WT = (df.loc[t, 24].values + df.loc[t, 25].values + df.loc[t, 26].values)/3

        for i in range(len(data_labels)):
            ave = [100*df.loc[t, 3*i+3].values/WT,
                   100*df.loc[t, 3*i+4].values/WT,
                   100*df.loc[t, 3*i+5].values/WT]
            flattened = np.hstack(ave).tolist()
            aves.append(flattened)

        """means = []
        stderrs = []

        for i in range(len(aves)):
            mean = np.average(aves[i])
            stderr = np.std(aves[i])/len(aves[i])
            means.append(mean)
            stderrs.append(stderr)

        df2 = pd.DataFrame({'mean':means, 'stderr':stderrs}, index=data_labels)
        df2 = df2.reset_index()
        df2.columns = ["label", "average", "stderr"]"""

        df2 = pd.DataFrame({"color" : colors, "protein" : proteins, "conc" : concentrations})
        #x = [df2, metadata_df]
        #df2 = pd.concat(x, axis=1)
        df2.conc = df2.conc.astype('category')
        df2.conc.cat.reorder_categories(natsorted(set(df2.conc)), inplace=True, ordered=True)
        df2.sort_values(by=["color", "protein", "conc"])

        lst1 = [item[0] for item in aves]
        lst2 = [item[1] for item in aves]
        lst3 = [item[2] for item in aves]

        df2["rep1"] = lst1
        df2["rep2"] = lst2
        df2["rep3"] = lst3

        """fig = plt.figure(figsize=figsize)
        plt.bar(range(len(df2)), df2.average, align='center', yerr=df2.stderr, color=colors, linewidth=1, edgecolor="black")
        plt.xticks(range(len(df2)), df2.label, rotation=60, fontsize=14)
        plt.yticks(fontsize=14)

        plt.plot([-0.5,int(len(df2)+0.5)], [100, 100], "k--")

        plt.ylabel("% growth at " + str(time_point) + " hours", fontsize=14)
        plt.title(plot_title, fontsize=16)
        plt.ylim(0, 125)
        plt.tight_layout()
        plt.show()"""

    if save == True:
        plt.savefig(output)

    return aves, df2
