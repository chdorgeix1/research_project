import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def VDBarChart(df):
    R1_list = []
    R2_list = []
    bar_chart_counter = 1
    for sample in df.samples.unique():
        df1 = df[df['samples'] == sample]
        #print(df1.head())
        df2 = df1[df1['paired_end'] == 'R1']
        #print(df2.head())
        for (columnName, columnData) in df2.iteritems():
            if columnName == 'samples' or columnName == 'paired_end':
                None
            else:
                R1_list.append(columnData.values[0])
        df3 = df1[df1['paired_end'] == 'R2']
        for (columnName, columnData) in df3.iteritems():
            if columnName == 'samples' or columnName == 'paired_end':
                None
            else:
                R2_list.append(columnData.values[0])

        a = 0 #is only in R1_list
        b = 0 #is in both lists
        c = 0 #is only in R2_list

        for i in range(len(R1_list)):   
            x = R1_list[i] - R2_list[i]
            if x == -1:
                c += 1
            elif x == 1:
                a += 1
            elif x == 0 and R1_list[i] == 1:
                b += 1

        #print(a, b, c)            
        #print(R1_list)
        #print(R2_list)
        
        # Plot the figure.

        if bar_chart_counter == 1:
            plt.bar(bar_chart_counter, a, color='b', label = 'R1 Only')
            plt.bar(bar_chart_counter, b, color='g', bottom=a, label = 'Both')
            plt.bar(bar_chart_counter, c, color='y', bottom=a + b, label = 'R2 Only')
        else:
            plt.bar(bar_chart_counter, a, color='b')
            plt.bar(bar_chart_counter, b, color='g', bottom=a)
            plt.bar(bar_chart_counter, c, color='y', bottom=a + b)
        # reordering the labels
        handles, labels = plt.gca().get_legend_handles_labels()

        # specify order
        order = [2, 1, 0]

        # pass handle & labels lists along with order as below
        plt.legend([handles[i] for i in order], [labels[i] for i in order])
        bar_chart_counter += 1
    ax = plt.axes()
    ax.set_xticks(list(range(1,len(df.samples.unique())+1)))
    ax.set_xticklabels((df.samples.unique()).tolist())

    plt.show()





