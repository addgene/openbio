import pandas as pd
import os, argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
​
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputfile',
        type=str,
        help='File with raw sequencing reads')
    parser.add_argument('--column',
        type=int,
        help='Column for which data has to be plotted')
    parser.add_argument('--outputfile',
        type=str,
        default='read_frac_plot.pdf',
        help='Path to outputfile')
    return parser
​
if __name__ == '__main__':
    args = get_parser().parse_args()
    inputfile = args.inputfile
    fig1,ax1 = plt.subplots()
    input_df = pd.read_table(inputfile, sep=',')
    c = args.column-1
    col_sum = float(sum(input_df.ix[:,c]))
    input_df['read_frac'] = [x/col_sum for x in input_df.ix[:,c]]
    input_df = input_df.sort_values(by='read_frac',ascending=False)
    input_df['Cumulative_sum'] = np.cumsum(input_df.read_frac)
    x_axis = [x/float(len(input_df)) for x in range(0,len(input_df))]
    y_axis = list(input_df.Cumulative_sum)
    ax1.plot(x_axis,y_axis,linewidth=1)
    ax1.set_xlim(0.0,1.0)
    ax1.set_ylim(0.0,1.0)
    auc = metrics.auc(x_axis,y_axis)
    ax1.text(0.6,0.2,'AUC = '+str(round(auc,2)),fontsize=14,fontweight='bold')
    ax1.tick_params(axis='both',labelsize=14,)
    ax1.set_xlabel('sgRNAs ranked by abundance',fontsize=14,fontweight='bold')
    ax1.set_ylabel('Cumulative fraction of total represented',fontsize=14,fontweight='bold')
    fig1.savefig(args.outputfile,format='pdf')
