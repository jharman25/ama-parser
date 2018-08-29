import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

def ama_parser_old_school(file_name, discard, blanks, data_labels, time_point, plot_title, sort_list,
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

    #rename time column in dataframe
    time = []

    for i in range(len(df)):
        if i == 0:
            time.append(i)
        else:
            time.append(time[i-1] + time_int)

    df = pd.DataFrame(d2).T

    df["Time (hours)"] = time

    #average triplicate data

    aves = []

    WT = (df.loc[:, 3].values + df.loc[:, 4].values + df.loc[:, 5].values)/3

    for i in range(len(data_labels)):
        ave = [100*df.loc[:, 3*i+3].values/WT,
               100*df.loc[:, 3*i+4].values/WT,
               100*df.loc[:, 3*i+5].values/WT]
        flattened = np.hstack(ave).tolist()
        aves.append(flattened)

    # calculate means

    means = []
    stderrs = []

    for i in range(len(aves)):
        mean = np.average(aves[i])
        stderr = np.std(aves[i])/len(aves[i])
        means.append(mean)
        means.append(mean)
        stderrs.append(stderr)

    df2 = pd.DataFrame({'mean':means, 'stderr':stderrs}, index=data_labels)
    df2 = df2.reset_index()
    df2.columns = ["protein", "average", "stderr"]
    df2['protein'] = pd.Categorical(df2['protein'], sort_list)
    df2 = df2.sort_values(by=['protein'])

    #optional plotting of triplicate-averaged growth curves when growth_curves == True
    if growth_curves == True:

        for i in range(len(data_labels)+1):

            means = []
            stddevs = []

            for j in range(len(df)):
                x = np.mean(df[df.columns[i*3+3:i*3+6]].values[j])
                y = np.std(df[df.columns[i*3+3:i*3+6]].values[j])
                means.append(x)
                stddevs.append(y)

            plt.errorbar(df["Time (hours)"], means, yerr=stddevs, fmt=".")
            if i == 0:
                plt.title("Wildtype growth")
            else:
                plt.title(data_labels[i-1])
            plt.xlabel("Time (hours)")
            plt.ylabel("OD600")
            plt.ylim(0,1.5)
            plt.show()

    #default boxplots of all data (leaves out blank)

    if plots == True:

        aves = []

        t = df.index[df["Time (hours)"] == time_point].values

        WT = (df.loc[t, 3].values + df.loc[t, 4].values + df.loc[t, 5].values)/3

        for i in range(len(data_labels)):
            ave = [100*df.loc[t, 3*i+6].values/WT,
                   100*df.loc[t, 3*i+7].values/WT,
                   100*df.loc[t, 3*i+8].values/WT]
            flattened = np.hstack(ave).tolist()
            aves.append(flattened)

        means = []
        stddevs = []

        for i in range(len(aves)):
            mean = np.average(aves[i])
            stddev = np.std(aves[i])
            means.append(mean)
            stddevs.append(stddev)

        df2 = pd.DataFrame({'mean':means, 'stddev':stddevs}, index=data_labels)
        df2 = df2.reset_index()
        df2.columns = ["protein", "average", "stddev"]
        df2['protein'] = pd.Categorical(df2['protein'], sort_list)
        df2 = df2.sort_values(by=['protein'])

        fig = plt.figure(figsize=figsize)
        plt.bar(range(len(df2)), df2.average, align='center', yerr=df2.stddev, color=colors, linewidth=1, edgecolor="black")
        plt.xticks(range(len(df2)), df2.protein, fontsize=14, rotation=60)
        plt.yticks(fontsize=14)

        plt.plot([-0.5,int(len(df2)+0.5)], [100, 100], "k--")

        plt.ylabel("% growth at " + str(time_point) + " hours", fontsize=14)
        plt.title(plot_title, fontsize=16)
        plt.ylim(0, 125)
        plt.tight_layout()
        plt.show()

    if save == True:
        plt.savefig(output)

    return df, df2
