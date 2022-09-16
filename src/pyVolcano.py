#!/usr/bin/env python
# coding: utf-8
"""Volcano plot from differential gene expression data.

This script allows the user to generate a volcano plot from a dataset in .csv,
.tsv or .xlsx format. The dataset must have a column for the gene names, the
p-values and the log2Fold changes. The figure is saved in the directory of
choice of the user with the given format (.pdf, .svg or any other format
admitted by matplotlib).

When used as a command line script, the following parameters should be given:

    Parameters
    ----------
    in_file : string.
        Path to the DE dataset. The dataset must contain the gene names, the p-values and the
        log2Fold values in columns.
    out_file : string.
        Path to the file where the plot is going to be saved.
    pval : float (optional).
        P-value threshold that determines significance. Defaults to 0.01.
    log2F : float (optional).
        Log2Fold value threshold that determines significance. Defaults to 1.
    gene_col : string (optional).
        Name of the column corresponding to gene names in the DE dataset. Defaults to 'gene'.
    pval_col : string (optional).
        Name of the column corresponding to p-values in the DE dataset. Defaults to 'padj'.
    log_col : string (optional).
        Name of the column corresponding to log2Fold values in the DE dataset. Defaults to
        'log2Fold'.
    n_names2show : int (optional).
        Number of top gene names to show. Defaults to 0.
    title : string (optional).
        Title of the plot to be written on top of the plot. Defaults to 'Volcano plot'.
    up_color : string (optional).
        Color for the up-regulated genes. Defaults to 'green'.
    down_color : string (optional).
        Color for the down-regulated genes. Defaults to 'red'.
    width : int (optional).
        Width of the figure in inches. Defaults to 8.
    height : int (optional).
        Height of the figure in inches. Defaults to 8.

The script must be run in the environment `containers/env.yml`, or following the
requirements in `requiremets.txt`.

This file can also be imported as a module and contains the following
functions:

    * volcano- returns the matplotlib.axes.Axes object with the plot.
    * main - the main function of the script.

Author: Juan Sebastian Diaz Boada
        juan.sebastian.diaz.boada@ki.se
07/07/22
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def volcano(df,ax,pval=0.01,log2F=1,gene_col='gene',pval_col='p-val',
            log_col='log2F',n_names2show = 10,title = 'Volcano plot',
                 up_color='green',down_color='red'):
    """Generates volcano plot from dataframe returning Axes with plot.

    Takes a differential gene expression dataset as a pandas.DataFrame (df) and
    generates a volcano plot highlighting over and under expressed genes, plus
    delimiting the significant areas of the plot. The plot is returned in a
    matplotlib.axes.Axes object.

    Parameters
    ----------
    df : pandas DataFrame.
        DataFrame holding the differential gene expression data with 3 columns:
        gene name, p-value and log2Fold change. Other columns are ignored but allowed.
    ax : matplotlib.axes.Axes.
        Axes where to plot the Volcano plot.
    pval : float (optional).
        P-value threshold that determines significance. Defaults to 0.01.
    log2F : float (optional).
        Log2Fold value threshold that determines significance. Defaults to 1.
    gene_col : string (optional).
        Name of the column corresponding to gene names in the DE dataset. Defaults to 'gene'.
    pval_col : string (optional).
        Name of the column corresponding to p-values in the DE dataset. Defaults to 'padj'.
    log_col : string (optional).
        Name of the column corresponding to log2Fold values in the DE dataset. Defaults to
        'log2Fold'.
    n_names2show : int (optional).
        Number of top gene names to show. Defaults to 0.
    title : string (optional).
        Title of the plot to be written on top of the plot. Defaults to 'Volcano plot'.
    up_color : string (optional).
        Color for the up-regulated genes. Defaults to 'green'.
    down_color : string (optional).
        Color for the down-regulated genes. Defaults to 'red'.

    Returns
    -------
    matplotlib.Axes.ax
        Ax object holding the plot.
    """
    DF = df.copy()
    # Parameter verification
    if pval<0 or pval>1:
        raise ValueError("P-value has to be between 0.0 and 1.0.")
    if log2F<0:
        raise ValueError("Log2Fold change has to be positive.")
    # Sort DF properly
    DF.insert(2,'absLogF',np.absolute(DF.loc[:,log_col]))
    DF = DF.sort_values([pval_col,log_col],ascending=[True,False]).reset_index(drop=True)
    # Insert color
    DF.insert(4,'color','black')
    down = (DF.loc[:,pval_col]<pval)&(DF.loc[:,log_col]<-log2F)
    up = (DF.loc[:,pval_col]<pval)&(DF.loc[:,log_col]>log2F)
    DF.loc[down,'color'] = down_color
    DF.loc[up,'color'] = up_color
    # Data variables for plot
    x = DF.loc[:,log_col].values
    y = -np.log10(DF.loc[:,pval_col].values)
    color = DF.loc[:,'color'].values
    # Plot
    ax.axhline(-np.log10(pval),color='gray',linestyle='--')
    ax.axvline(log2F,color='gray',linestyle='--')
    ax.axvline(-log2F,color='gray',linestyle='--')
    ax.scatter(x,y,c=color,s=3)
    # Truncation of names to show
    DE = DF.loc[(DF['absLogF']>log2F) & (DF[pval_col]<pval),:]
    n_sg = len(DE)
    n_names2show = n_names2show if n_names2show<n_sg else n_sg
    names = DE.loc[:,gene_col].values
    x_de = DE.loc[:,log_col].values
    y_de = -np.log10(DE.loc[:,pval_col].values)
    for i in range(n_names2show):
        ha = "right" if x_de[i] > 0 else "left"
        ax.text(x_de[i], y_de[i] , s = names[i],ha=ha)
    ax.set_ylabel(r'$-log_{10}(pval)$')
    ax.set_xlabel(r'$log_2FoldChange$')
    ax.set_title(title)
    return ax

def main():
    # Argparser definition
    import argparse
    parser = argparse.ArgumentParser(description="Parameters of volcano plot.")
    parser.add_argument('in_file', type=str, help="Path to the DE dataset.")
    parser.add_argument('out_file', type=str, help="Path to the file where the figure will be saved.")
    parser.add_argument('--pval', type=float, default=0.01, help="P-value threshold to determine significance. Defaults to 0.01.")
    parser.add_argument('--log2F', type=float, default=1, help="Log2Fold threshold to determine significance. Defaults to 1.")
    parser.add_argument('--gene_col', type=str, default='gene', help="Name of the column corresponding to gene names in the DE dataset. Defaults to 'gene'.")
    parser.add_argument('--pval_col', type=str, default='p-val', help="Name of the column corresponding to p-values in the DE dataset. Defaults to 'padj'.")
    parser.add_argument('--log_col', type=str, default='log2Fold', help="Name of the column corresponding to log2 values in the DE dataset. Defaults to 'log2FoldChange'.")
    parser.add_argument('-n','--n_names2show', type=int, default=0, help="Number of top gene names to show. Defaults to 0.")
    parser.add_argument('--title', type=str, default='Volcano plot', help="Title of the plot. Defaults to 'Volcano plot'.")
    parser.add_argument('--up_color', type=str, default='green', help="Color for up-regulated genes. Defaults to 'green'.")
    parser.add_argument('--down_color', type=str, default='red', help="Color for down-regulated genes. Defaults to 'red'.")
    parser.add_argument('--width', type=int, default=8, help="Width of the figure in inches. Defaults to 8.")
    parser.add_argument('--height', type=int, default=8, help="Height of the figure in inches. Defaults to 8.")
    args = parser.parse_args()
    # Assignment of parameter variables
    in_file = args.in_file
    out_file = args.out_file
    pval = args.pval
    log2F = args.log2F
    gene_col = args.gene_col
    pval_col = args.pval_col
    log_col = args.log_col
    n_names2show = args.n_names2show
    title = args.title
    up_color = args.up_color
    down_color = args.down_color
    width = args.width
    height = args.height

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
    fig,ax = plt.subplots(1,1,figsize=[width,height])
    ax = volcano(DF,ax,pval=pval,log2F=log2F,gene_col=gene_col,pval_col=pval_col,
                 log_col=log_col,n_names2show = n_names2show,title = title,
                 up_color=up_color,down_color=down_color)
    plt.savefig(out_file)

if __name__ == "__main__":
    main()
