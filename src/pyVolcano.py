#!/usr/bin/env python
# coding: utf-8

# # pyVolcano
# Volcano plot working over matplotlib, numpy and pandas

import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# Argparser definition
parser = argparse.ArgumentParser(description="Parameters of volcano plot.")
parser.add_argument('in_path', type=str, help="Path to the DE dataset.")
parser.add_argument('out_path', type=str, help="Path to the file where the figure will be saved.")
parser.add_argument('--pval', type=float, default=0.01, help="P-value threshold to determine significance. Defaults to 0.01.")
parser.add_argument('--log2F', type=float, default=1, help="Log2Fold threshold to determine significance. Defaults to 1.")
parser.add_argument('-n','--n_names2show', type=int, default=0, help="Number of top gene names to show. Defaults to 0.")
parser.add_argument('--pval_col', type=str, default='padj', help="Name of the column corresponding to p-values in the DE dataset. Defaults to 'padj'.")
parser.add_argument('--log_col', type=str, default='log2FoldChange', help="Name of the column corresponding to log2 values in the DE dataset. Defaults to 'padj'.")
parser.add_argument('--gene_col', type=str, default='gene', help="Name of the column corresponding to gene names in the DE dataset. Defaults to 'padj'.")
parser.add_argument('--title', type=str, default='Volcano plot', help="Title of the plot. Defaults to 'Volcano plot'.")
parser.add_argument('--up_color', type=str, default='green', help="Color for up-regulated genes. Defaults to 'green'.")
parser.add_argument('--down_color', type=str, default='red', help="Color for down-regulated genes. Defaults to 'red'.")
args = parser.parse_args()
# Assignment of parameter variables
in_file = args.in_path
out_file = args.out_path
pval_thresh = args.pval
log_thresh = args.log2F
pval_col = args.pval_col
log_col = args.log_col
gene_col = args.gene_col
title = args.title
n_names2show = args.n_names2show
up_color = args.up_color
down_color = args.down_color

# Import dataset
DF = pd.read_csv(in_file,sep='\t',index_col=0).reset_index(drop=True)

# Sort DF properly

DF.insert(2,'absLogF',np.absolute(DF.loc[:,'log2FoldChange']))
DF = DF.sort_values(['padj','log2FoldChange'],ascending=[True,False]).reset_index(drop=True)

# Insert color
DF.insert(4,'color','black')
down = (DF.loc[:,pval_col]<pval_thresh)&(DF.loc[:,log_col]<-log_thresh)
up = (DF.loc[:,pval_col]<pval_thresh)&(DF.loc[:,log_col]>log_thresh)
DF.loc[down,'color'] = down_color
DF.loc[up,'color'] = up_color

# Plot

x = DF.loc[:,log_col].values
y = -np.log10(DF.loc[:,pval_col].values)
names = DF.loc[:,gene_col]

fig,ax = plt.subplots(figsize=[8,8])
ax.axhline(-np.log10(pval_thresh),color='gray',linestyle='--')
ax.axvline(log_thresh,color='gray',linestyle='--')
ax.axvline(-log_thresh,color='gray',linestyle='--')
ax.scatter(x,y,c=DF.loc[:,'color'].values,s=3)
for i in range(n_names2show):
    ha = "right" if x[i] > 0 else "left"
    ax.text(x[i], y[i] , s = names[i],ha=ha)
ax.set_ylabel(r'$-log_{10}(pval)$')
ax.set_xlabel(r'$log_2FoldChange$')
ax.set_title(title)
plt.savefig(out_file,format='pdf')

# In[ ]:
