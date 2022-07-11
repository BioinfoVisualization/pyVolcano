#!/usr/bin/env python
# coding: utf-8
# =========================================================================== #
# pyVolcano.py                                                                #
# Author: Juan Sebastian Diaz Boada                                           #
# Creation Date: 07/07/22                                                     #
# =========================================================================== #

""" Creates a volcano plot from a dataset and exports it.
    Reads a dataset comming from a differential expression experiment, containing a list of genes
    with their respective p-values and log2Fold values and generates a volcano plot. The plot is
    exported in a file and format given by the user.

    Parameters
    ----------
    in_path : string.
        Path to the DE dataset. The dataset must contain the gene names, the p-values and the
        log2Fold values in columns.
    out_path : string.
        Path to the file where the plot is going to be saved.
    pval : float (optional).
        P-value threshold that determines significance. Defaults to 0.01.
    log2F : float (optional).
        Log2Fold value threshold that determines significance. Defaults to 1.
    n_names2show : int (optional).
        Number of top gene names to show. Defaults to 0.
    pval_col : string (optional).
        Name of the column corresponding to p-values in the DE dataset. Defaults to 'padj'.
    log_col : string (optional).
        Name of the column corresponding to log2Fold values in the DE dataset. Defaults to
        'log2FoldChange'.
    gene_col : string (optional).
        Name of the column corresponding to gene names in the DE dataset. Defaults to 'gene'.
    title : string (optional).
        Title of the plot to be written on top of the plot. Defaults to 'Volcano plot'.
    up_color : string (optional).
        Color for the up-regulated genes. Defaults to 'green'.
    down_color : string (optional).
        Color for the down-regulated genes. Defaults to 'red'.

"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
#---------------------------------------------------------------------------------------------------#
def volcano(df,ax,gene_col='gene',pval_col='p-val',log_col='log2F',n_names2show = 10,
                 pval_thresh=0.01,log_thresh=1,title = 'Volcano plot',
                 up_color='green',down_color='red'):
    DF = df.copy()
    # Sort DF properly
    DF.insert(2,'absLogF',np.absolute(DF.loc[:,log_col]))
    DF = DF.sort_values([pval_col,log_col],ascending=[True,False]).reset_index(drop=True)
    # Truncation of names to show
    n_sg = len(DF.loc[(DF['absLogF']>log_thresh) & (DF[pval_col]<pval_thresh),:])
    n_names2show = n_names2show if n_names2show<n_sg else n_sg
    # Insert color
    DF.insert(4,'color','black')
    down = (DF.loc[:,pval_col]<pval_thresh)&(DF.loc[:,log_col]<-log_thresh)
    up = (DF.loc[:,pval_col]<pval_thresh)&(DF.loc[:,log_col]>log_thresh)
    DF.loc[down,'color'] = down_color
    DF.loc[up,'color'] = up_color
    # Data variables for plot
    x = DF.loc[:,log_col].values
    y = -np.log10(DF.loc[:,pval_col].values)
    color = DF.loc[:,'color'].values
    names = DF.loc[:,gene_col]
    # Plot
    ax.axhline(-np.log10(pval_thresh),color='gray',linestyle='--')
    ax.axvline(log_thresh,color='gray',linestyle='--')
    ax.axvline(-log_thresh,color='gray',linestyle='--')
    ax.scatter(x,y,c=color,s=3)
    for i in range(n_names2show):
        ha = "right" if x[i] > 0 else "left"
        ax.text(x[i], y[i] , s = names[i],ha=ha)
    ax.set_ylabel(r'$-log_{10}(pval)$')
    ax.set_xlabel(r'$log_2FoldChange$')
    ax.set_title(title)
    return ax
#---------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    # Argparser definition
    import argparse
    parser = argparse.ArgumentParser(description="Parameters of volcano plot.")
    parser.add_argument('in_path', type=str, help="Path to the DE dataset.")
    parser.add_argument('out_path', type=str, help="Path to the file where the figure will be saved.")
    parser.add_argument('--pval', type=float, default=0.01, help="P-value threshold to determine significance. Defaults to 0.01.")
    parser.add_argument('--log2F', type=float, default=1, help="Log2Fold threshold to determine significance. Defaults to 1.")
    parser.add_argument('-n','--n_names2show', type=int, default=0, help="Number of top gene names to show. Defaults to 0.")
    parser.add_argument('--gene_col', type=str, default='gene', help="Name of the column corresponding to gene names in the DE dataset. Defaults to 'gene'.")
    parser.add_argument('--pval_col', type=str, default='p-val', help="Name of the column corresponding to p-values in the DE dataset. Defaults to 'padj'.")
    parser.add_argument('--log_col', type=str, default='log2Fold', help="Name of the column corresponding to log2 values in the DE dataset. Defaults to 'log2FoldChange'.")
    parser.add_argument('--title', type=str, default='Volcano plot', help="Title of the plot. Defaults to 'Volcano plot'.")
    parser.add_argument('--up_color', type=str, default='green', help="Color for up-regulated genes. Defaults to 'green'.")
    parser.add_argument('--down_color', type=str, default='red', help="Color for down-regulated genes. Defaults to 'red'.")
    args = parser.parse_args()
    # Assignment of parameter variables
    in_file = args.in_path
    out_file = args.out_path
    pval_thresh = args.pval
    log_thresh = args.log2F
    n_names2show = args.n_names2show
    gene_col = args.gene_col
    pval_col = args.pval_col
    log_col = args.log_col
    title = args.title
    up_color = args.up_color
    down_color = args.down_color

    # Import dataset
    file_type = in_file.split('.')[-1]
    if file_type == 'tsv':
        DF = pd.read_csv(in_file,sep='\t',index_col=0).reset_index(drop=True)
    elif file_type == 'xlsx':
        DF = pd.read_excel(in_file,index_col=0)
    elif file_type == 'csv':
        DF = pd.read_csv(in_file,sep=',',index_col=0).reset_index(drop=True)
    else:
        raise NameError("Invalid input format. Has to be either .tsv, .csv or .xlsx.")

    # Plot
    fig,ax = plt.subplots(1,1,figsize=[8,8])
    ax = volcano(DF,ax,gene_col=gene_col,pval_col=pval_col,log_col=log_col,
                    n_names2show = n_names2show,pval_thresh=pval_thresh,
                    log_thresh=log_thresh,title = title,up_color=up_color,
                    down_color=down_color)
    plt.savefig(out_file)
