# pyVolcano
Volcano plot in python! :volcano:<br>
`pyVolcano` is a small module build over `numpy`, `matplotlib` and `pandas` that creates a Volcano plot from a dataset using the function [`volcano`](https://github.com/BioinfoVisualization/pyVolcano/blob/5770249d04e28d9d2e96281da754753ce5b38ad1/src/pyVolcano.py#L61)

## 1. Installation
`pyVolcano` can be used in 3 different ways.
### 1.1 Installation of packages in [requirements.txt](requirements.txt)
`pyVolcano` is tested in Python 3.8 and its dependencies can be pip-installed:
```bash
pip install -r requirements.txt
```
### 1.2 Setting a conda environment
This repository contains a conda definition file [`env.yml`](containers/env.yml) that can create the conda environment:
```bash
conda env create -f env.yml
conda activate volcano
```
### 1.3 Singularity container
For Linux users, `pyVolcano` can be run easily with a Singularity container. The container can be built with the definition file [`pyVolcano.def`](containers/pyVolcano.def) and run as follows:
```bash
sudo singularity build pyVolcano.sif containers/pyVolcano.def
./pyVolcano.sif <module_parameters>
```
if the user does not have admin rights, the `--fakeroot` flag can be used when building the singularity image:
```bash
singularity build --fakeroot pyVolcano.sif containers/pyVolcano.def
./pyVolcano.sif <module_parameters>
```
## 2. Usage
`pyVolcano` can be used as a command-line tool running [`pyVolcano.py`](pyVolcano.py) or by loading it as a python module and using its function [`volcano`](https://github.com/BioinfoVisualization/pyVolcano/blob/5770249d04e28d9d2e96281da754753ce5b38ad1/src/pyVolcano.py#L61).
### 2.1 Using the command line
By calling the script with the corresponding parameters in the correct environment, the script will save the plot in the location given by the user and in the given format.
```bash
./pyVolcano.py in_file out_file --pval 0.01 --log2F 1 --gene_col gene
               --pval_col p-value --log_col log2F -n 0 --title 'Volcano plot'
               --up_color green --down_color red --width 8 --height 8
```
#### Parameters for command-line calling
+ `in_file` : string. Path to the DE dataset. The dataset must contain the gene names, the p-values and the log2Fold values in columns as follows:<br>

|gene | log2Fold | p-value|
|------------ | ------------- | -------------|
|ENSG00000162391 |-2.569 | 0.0054|
|ENSG00000181965 | 1.559 | 0.0015|
... | ... | ...|

The names in the column *genes* are the ones to be shown in the plot for selected genes (see below `n_names2show`). The name of the columns do not need to be the same as in this example (see below `gene_col`, `pval_col` and `log_col`).
+ `out_file` : string. Path to the file where the plot is going to be saved. The name of the file must have the desiderd format. Available formats are the ones supported by `matplotlib`, such as `.pdf`, `.png`, `.svg`, etc.
+ `pval` : float (optional). P-value threshold that determines significance. A horizontal line will be drawn in the plot corresponding to this value. Defaults to 0.01.
+ `log2F` : float (optional). Log2Fold value threshold that determines significance. Two vertical lines will be drawn in the plot corresponding to this value and its opposite. Defaults to 1.
+ `gene_col` : string (optional). Name of the column corresponding to gene names in the DE dataset. Defaults to 'gene'.
+ `pval_col` : string (optional). Name of the column corresponding to p-values in the DE dataset. Defaults to 'padj'.
+ `log_col` : string (optional). Name of the column corresponding to log2Fold values in the DE dataset. Defaults to 'log2Fold'.
+ `n_names2show` : int (optional). Number of top gene names to show. The genes are sorted by ascending p-value and by descending absolute value of log2Fold change, giving priority to the former. If the number of significant genes is higher than this number, only the names of significant genes are shown. Defaults to 0.
+ `title` : string (optional). Title of the plot to be written on top of the plot. Defaults to 'Volcano plot'.
+ `up_color` : string (optional). Color for the up-regulated genes. Defaults to 'green'.
+ `down_color` : string (optional). Color for the down-regulated genes. Defaults to 'red'.
+ `width` : int (optional). Width of the figure in inches. Defaults to 8.
+ `height` : int (optional). Height of the figure in inches. Defaults to 8.
### 2.2 Importing as a module
A script or notebook running in the same directory as [`pyVolcano.py`](pyVolcano.py) can import the function [`volcano`](https://github.com/BioinfoVisualization/pyVolcano/blob/5770249d04e28d9d2e96281da754753ce5b38ad1/src/pyVolcano.py#L61) as follows:
```python
from pyVolcano import volcano
```
As the function [`volcano`](https://github.com/BioinfoVisualization/pyVolcano/blob/5770249d04e28d9d2e96281da754753ce5b38ad1/src/pyVolcano.py#L61) returns a `matplotlib.axes.Axes` object, a figure and an axes should be created to make the plot.
```python
from matplotlib import pyplot as plot
fig,ax = plt.subplots()
ax = volcano(df,ax,pval=0.01,log2F=1,gene_col='gene',pval_col='p-val',
        log_col='log2F',n_names2show = 10,title = 'Volcano plot',
        up_color='green',down_color='red')
plt.show()
```
#### Parameters of python function
+ `df` : pandas DataFrame holding the differential gene expression data with the same structure as the input file explained above. Other columns are ignored but allowed.
+ `ax` : matplotlib.axes.Axes. Axes where to plot the Volcano plot.
+ `pval` : float (optional). P-value threshold that determines significance. Defaults to 0.01.
+ `log2F` : float (optional). Log2Fold value threshold that determines significance. Defaults to 1.
+ `gene_col` : string (optional). Name of the column corresponding to gene names in the DE dataset. Defaults to 'gene'.
+ `pval_col` : string (optional). Name of the column corresponding to p-values in the DE dataset. Defaults to 'padj'.
+ `log_col` : string (optional). Name of the column corresponding to log2Fold values in the DE dataset. Defaults to 'log2Fold'.
+ `n_names2show` : int (optional).Number of top gene names to show. Defaults to 0.
+ `title` : string (optional). Title of the plot to be written on top of the plot. Defaults to 'Volcano plot'.
+ `up_color` : string (optional). Color for the up-regulated genes. Defaults to 'green'.
+ `down_color` : string (optional). Color for the down-regulated genes. Defaults to 'red'.
