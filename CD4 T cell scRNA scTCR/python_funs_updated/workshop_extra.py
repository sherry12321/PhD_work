import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import anndata
import seaborn as sns
from anndata import AnnData
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import coo_matrix, csr_matrix
from typing import Optional, Union, Mapping  # Special
from typing import Sequence  # ABCs
from typing import Tuple  # Classes

def compute_entropy(adata, obs_index: str = None, neighbors_index: str = 'neighbors', k: int = 30, pca_key: str = None, inplace=True, key_added: str = 'entropy'):
    
    '''
      Entropy analysis:
      adata: Anndata object
      obs_index: key for sample IDs categories (ie 'Timepoints' or 'batch_id')
      neighbors_index: place to look for nearest neighbors if already computed
      k = # nearest neighbors to compute if not yet computed
      pca_key = pca key to run neighbor computation on if not computed already
      inplace - whether to add the results to obs in adata or return the pd dataframe
      key_added - for where to store the entropy analysis info (if inplace = True)
  '''
    # Preprocessing checks:
    if neighbors_index not in adata.uns:  # if needed, run neighbors
        sc.pp.neighbors(adata, n_neighbors = k, use_rep = pca_key)
        neighbor_arr = adata.obsp['distances']  # sparse matrix of distances
    elif neighbors_index == 'neighbors':
        neighbor_arr = adata.obsp['distances']
    else:
        neighbor_arr = adata.obsp[neighbors_index + '_distances']

    if obs_index not in adata.obs.columns:
        raise ValueError('No batch ids provided and no obs key index provided')
    else:
        sample_ids = adata.obs[obs_index]  # contains batch ids

    # Build n x k arr of cells x k-neighbor indices
    n_neighbors = adata.uns[neighbors_index]['params']['n_neighbors']
    num_cells = adata.shape[0]
    # n x k-1 neighbors-doesn't include self
    neighbor_inds = np.reshape(neighbor_arr.indices, newshape=(num_cells, n_neighbors-1))
    cell_ints = np.arange(num_cells, dtype=int)  # add self to matrix of k nearest neighbors
    neighbor_inds = np.column_stack((cell_ints, neighbor_inds))

    # Map batch ids to ints
    batches = sample_ids.unique()
    batch_map = {}  # dict of batch, int representing id in order of batches
    for batch_ind in range(len(batches)):
        batch_map[batches[batch_ind]] = batch_ind
    cell_to_batch = np.array([batch_map[sample_id] for sample_id in sample_ids])  # n cells -> batch id int

    # Compute entropy
    output_df = pd.Series(index = adata.obs.index, dtype='float64')
    for i in range(num_cells):
        # replace cell index w batch index, get unique+counts
        _, cts = np.unique(cell_to_batch[neighbor_inds[i]], return_counts=True)
        batch_freq = cts/n_neighbors
        cell_entropy = -(sum((batch_freq)*np.log2(batch_freq)))  # Shannon's Entropy
        output_df[i] = cell_entropy

    # Save as part of anndata if inplace
    if inplace:
        key_k = key_added
        adata.obs[key_k] = output_df

    if not inplace:
        return output_df, neighbor_inds



def volcano_plot(
    diff_expr_result : pd.DataFrame,
    logfoldchange_label: str,
    pvalue_label: str,
    logfold_neg_thres: Optional[float] = None,
    logfold_pos_thres: Optional[float] = None,
    pval_thres: Optional[float] = None,
    c: Optional[str] = 'lightgray',
    pos_c: Optional[str] = 'b',
    neg_c: Optional[str] = 'r',
    s: Optional[float] = 4,
    genes_highlight: Optional[list] = None,
    genes_highlight_fontdict: Optional[dict] = {'fontsize': 16},
    fig_size: Optional[tuple] = (8, 6),
    savefig: Optional[bool] = False,
    saveas: Optional[str] = "volcanoplot.png",
):
    """
    Plot ranking of genes using volano plot

    Visualizes and identifies statistifically significant gene expression changes from two different
    experimental conditions in terms of log fold change and p value

    Parameters
    ----------
    diff_expr_result
        pandas dataframe containing logfoldchanges, pvalues.
        The index must be gene names if want to highlight specific genes (see genes_highlight parameter)
    logfoldchange_label
        column name for log fold change in diff_expr_result
    pvalue_label
        column name for pvalues in diff_expr_result
    logfold_neg_thres, logfold_pos_thres and pval_thres
        Threshold for logfoldchange and pvalues to identify significant genes
    c
        color for the dots
    pos_c
        color for dots representing the upregulated genes
    neg_c
        color for dots representing the downregulated genes
    s
        dot size
    genes_highlight
        genes to highlight (must be a list of less than 10 genes)
    genes_highlight_fontdict
        additional properties for the gene names to be provided as a dictionary
    fig_size
        size of the figure, a tuple
    savefig
        boolean to save figure or not
    saveas
        Full path and name of output plot image (include ".png" extension)

    Returns
    -------
    Returns axes and saves figure
    
      
    Example
    -------
    >> f = volcano_plot(diff_expr_result = df,
                        logfoldchange_label = 'coef',
                        pvalue_label = 'fdr',
                        logfold_neg_thres = -0.3,
                        logfold_pos_thres = 0.3,
                        pval_thres = 10,
                        s = 10,
                        genes_highlight = ['TTR', 'PENK', 'TNF'],
                        genes_highlight_fontdict = {'fontsize': 20, 'fontstyle': 'italic'}, 
                        fig_size = (8, 6))
    f.legend(['ns', 'u', 'd'], markerscale = 2)
    f.axes.set_xlabel('log fold change', fontsize = 16)
    f.axes.set_ylabel('-log10(pvalues)', fontsize = 16)
    f.set_xticklabels(f.get_xticklabels(), fontsize = 16);
    f.set_yticklabels(f.get_yticklabels(), fontsize = 16);
    f.set_xlim([-0.85, 0.85])
    
    """
    
    logfoldchange = diff_expr_result[logfoldchange_label].values
    p_vals = diff_expr_result[pvalue_label].values
    
    # What to do with 0 p-value: we reassign it as minimum non-zero value divided by 2
    p_vals[p_vals == 0] = np.min(p_vals[p_vals != 0])/2
    logpvals = np.log10(p_vals)
    #logpvals[logpvals == -np.inf] = np.min(logpvals[logpvals != -np.inf])
    logpvals = -1*logpvals
    
    fig = plt.figure(figsize = fig_size)
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(logfoldchange, logpvals, c = c, s = s, label = 'Not Significant')
    if pval_thres:
        # find genes with significant expression
        sig_cond_pos = (logpvals > pval_thres) & (logfoldchange > logfold_pos_thres)
        sig_cond_neg = (logpvals > pval_thres) & (logfoldchange < logfold_neg_thres)
        id_sig_pos = np.where(sig_cond_pos)[0]
        id_sig_neg = np.where(sig_cond_neg)[0]
        
        ax.scatter(logfoldchange[id_sig_pos], logpvals[id_sig_pos], c = pos_c, s = s, label = 'Up-regulated')
        ax.scatter(logfoldchange[id_sig_neg], logpvals[id_sig_neg], c = neg_c, s = s, label = 'Down-regulated')
        
    if genes_highlight:
        if len(genes_highlight) < 10:
            genes_id = [diff_expr_result.index.get_loc(j) for j in genes_highlight]
            
            for k_iter, item in enumerate(genes_id):
                ax.text(logfoldchange[item], logpvals[item], genes_highlight[k_iter],
                    fontdict = genes_highlight_fontdict)
    
    ax.legend()
    if savefig:
        fig.savefig(saveas, bbox_inches = 'tight')
    return ax
    
    
def compute_mnn(arr1: np.ndarray, arr2: np.ndarray, k1: int = 30, k2: int = 30, dist_metric: str = 'euclidean', connectivities: bool = False):
    """ Computes the mutual nearest neighbors between two datasets
    Takes in two subsets of data or datasets across which to compute
    mutual nearest neighbors. Can compute different numbers of neighbors
    in each direction. Returns connectivities sparse matrix and distances
    sparse matrix, not guaranteed to have any k non-zero values for any
    row/column.
    Parameters
    ----------
    arr1
      First subset of the data to compute mutual nearest neighbors of (often
      in PCA space). For n1 number of cells, and k dimensions, this will
      have shape n1 x k
    arr2
      Second subset of the data for the mutual nearest neighbor computation,
      which should be in the same dimensionality space as arr1 (whether PCA
      or some other space). For n2 number of cells, and the same k dimensions, this
      will have shape n2 x k
    k1
      For each cell in arr2, find the k1 nearest neighbors from the arr1 data.
      Default value 30.
    k2
      For each cell in arr1, find the k2 nearest neighbors from the arr2 data.
      Default value is 30.
    dist_metric
      Distance metric to compute neighbors, from sklearn. Default is ``euclidean``.
      Valid metrics are: ``euclidean``, ``manhattan``, ``chebyshev``, ``minkowski``,
      ``wminkowski``, ``seuclidean``, ``mahalanobis``, ``hamming``, ``canberra``,
      ``braycurtis``
    connectivities
      Whether to return the connectivities matrix, otherwise distance matrix is returned. Default `False`

    Returns
    -------
    :class: `scipy.sparse.csr.csr_matrix`
      If connectivities, a matrix of dimension n1 x n2 with 1 in cells (i,j) if n1_i and
      n2_j are mutual nearest neighbors, 0 elsewhere.
      If connectivities is false, a matrix of dimension n1 x n2 where cell (i,j) is the
      distance between n1_i and n2_j if they are mutual nearest neighbors, otherwise 0.
  """
    # sklearn NearestNeighbors objects for each input array
    nbrs1 = NearestNeighbors(n_neighbors=k1, metric=dist_metric)
    nbrs2 = NearestNeighbors(n_neighbors=k2, metric=dist_metric)

    # Fit nbrs1 object to arr1, query k1 neighbors from arr2
    nbrs1.fit(arr1)
    t1_nbrs = nbrs1.kneighbors_graph(arr2, mode='distance')

    # Fit nbrs2 object to arr2, query k2 neighbors from arr1
    nbrs2.fit(arr2)
    t2_nbrs = nbrs2.kneighbors_graph(arr1, mode='distance')

    # Mututally nearest neighbors
    mnn = t2_nbrs.multiply(t1_nbrs.T)  # anywhere not mutual will have value 0
    mnn = mnn.sqrt()  # correct distances post-multiplication

    if connectivities:  # return connectivity matrix instead of distance matrix
        n_rows, n_cols = mnn.shape
        # calculate connectivities:
        row_connectivity, col_connectivity = mnn.nonzero()
        data = np.ones(len(row_connectivity))
        connectivities_matrix = csr_matrix((data, (row_connectivity, col_connectivity)), shape=(n_rows, n_cols))
        return connectivities_matrix

    return mnn


def check_NA_MAST(diff_expr_file, adata, groupby, group1, group2 = 'rest'):
    """ Sometimes you will see NA in log-fold change (coef) for specific genes in your differential expression results using MAST. NA will occur usually when the specific gene has 0 expression in either of the two groups you are comparing. For example, if you are computing differential expression between cluster-0 vs rest of the data and observe that Gene-X has NA in coef but has significant p-value then it is likely that Gene-X has super high expression in Cluster-0 but no expression in rest of the data. This is definitely a differentially expressed genes (hence the significant p-value) but MAST cannot compute log-fold change because there is 0 expression in the rest of the data, which leads to the NA. To answer this question, we need to know what the use of the differential analysis is for. (1) If you want to use the differential analysis results to annotate a cluster then this may not make much difference because this takes place rarely. (2) If you want to use the differential analysis results to perform Gene Set Enrichment Analysis then all it matters is the sign of coef and not the magnitude. This is a simple function to check for the cause of NA and get at least the sign of the coef.
    
    Parameters
    ----------
    diff_expr_file
      Path to a csv file of differential expression result
    adata
      AnnData object on which differential expression analysis was performed. It will be used to investigate expression of genes
    groupby
      Name of the clustering result (typically in adata.obs) that was used for differential expression analysis
    group1
      Identifier (e.g. cluster name) for the first group in differential analysis performed
    group2
        Identifier (e.g. cluster name) for the first group in differential analysis performed. By default it is set to "rest" indicating all the cells in the data except those in group1.
        
    Returns
    -------
    Returns edited MAST result (with no NA) and the dataframe (indexed by genes with NA in the original results) with their expression in each group.
    """

    # load the differential expression analysis result
    mast_result = pd.read_csv(diff_expr_file, index_col = 0)
    
    # identify which cells belong to group 1
    cells1 = np.where(adata.obs[groupby] == group1)[0]
    
    # identify which cells belong to group 2
    if group2 == 'rest':
        cells2 = np.where(adata.obs[groupby] != group1)[0]
    else:
        cells2 = np.where(adata.obs[groupby] == group2)[0]
    
    # identify the genes with NA in the differential expression result
    genes_with_NA = mast_result['primerid'][np.isnan(mast_result['coef'])].values
    
    # identify the location of the genes with NA
    genes_with_NA_index = [adata.var_names.get_loc(j) for j in genes_with_NA]
    
    # compute median expression of the genes in each of the groups
    expr_1 = np.asarray(np.mean(adata.X[cells1, :][:, genes_with_NA_index], axis = 0)).flatten()
    expr_2 = np.asarray(np.mean(adata.X[cells2, :][:, genes_with_NA_index], axis = 0)).flatten()
    
    # compute the sign of the difference in median expression
    df_temp = pd.DataFrame({'expr_1': expr_1, 'expr_2': expr_2, 'diff': expr_1 - expr_2}, index = list(genes_with_NA))
    df_temp['sign'] = [np.sign(j) for j in df_temp['diff']]

    # assign the sign back into mast results
    for gene_name in df_temp.index:
        id_loc = mast_result.index[mast_result['primerid'] == gene_name]
        mast_result.loc[id_loc, 'coef'] = df_temp.loc[gene_name]['sign']
        
    return mast_result, df_temp
    
    
    
    

def gsea_dotplot(df: pd.DataFrame,
                 xlab: str,
                 ylab: str,
                 hue = 'k',
                 log10_pvalue: Optional[bool] = False,
                 fig_size: Optional[tuple] = (12, 6),
                 line_color: Optional[str] = 'skyblue',
                 dot_size: Optional[float] = 100,
                 cmap: Optional[str] = 'viridis',
                 grid_on: Optional[bool] = True,
                 label_fontsize: Optional[float] = 16,
                 tick_fontsize: Optional[float] = 16,
                 colorbar: Optional[bool] = False,
                 savefig: Optional[Union[bool, str]] = False):
    
    
    """
    Make a dotplot to visualize GSEA results

    Visualizes specific pathways and associated parameters (expected output of GSEA)

    Parameters
    ----------
    df
        pandas dataframe containing pathway names, NES values and p-values.
    xlab
        column name in df for NES result (on x-axis)
    ylab
        column name in df for pathway or term (on y-axis)
    hue
        Either a color name or column name in df for coloring the dots (typically FDR q-val)
    log10_pvalue
        Boolean variable indicating whether to log10 hue (pvalue or fdr pvalue) or not. Default: False
    fig_size
        A tuple indicating size of the figure. Default: (12, 6)
    line_color
        Color of the horizontal lines. Default: 'skyblue'
    dot_size
        Size of the dots at the end of the horizontal lines. Default: 100
    cmap
        If hue indicates a continuous variable then choose colormap. Default: viridis
    grid_on
        Boolean indicating if grid lines are required in the plot. Default: True
    label_fontsize
        Size of the x-axis, y-axis and colorbar labels. Default: 16
    tick_fontsize
        Size of the x-axis, y-axis and colorbar tick labels. Default: 16
    colorbar
        Boolean to make a colorbar or not. Default: False
    savefig
        If you want to save figure, provide full path and name of output plot image (include ".png" extension). Default: False

    Returns
    -------
    Returns axes and saves figure
    
      
    Example
    -------
    gsea_dotplot(df,
                 xlab = 'NES',
                 ylab = 'Term',
                 hue = 'FDR q-val',
                 log10_pvalue = False,
                 colorbar = True,
                 tick_fontsize = 16,
                 label_fontsize = 16,
                 savefig = '/Users/sharmar1/Dropbox/msk_workshop/materials_2023_2024_batch1/session-5/gseapy_dotplot.png')
    
    """
    
    # check if hue exists in df
    if hue in df.columns:
        colorby = df[hue].values
        colorbar_label = hue
        if log10_pvalue:
            colorby[colorby == 0] = np.min(colorby[colorby != 0])/2
            colorby = -np.log10(colorby.astype(float))
            colorbar_label = '-log10(' + hue + ')'
    else:
        colorby = hue
        cmap = None
    

    # y-axis
    my_range = np.arange(1,len(df.index)+1)

    # create a figure
    fig = plt.figure(figsize = fig_size)
    ax = fig.add_subplot(1, 1, 1)

    # add horizontal lines
    ax.hlines(y=my_range, xmin=0, xmax=df[xlab], color=line_color)

    # add the circles
    im1 = ax.scatter(df[xlab], my_range, s = dot_size, c = colorby, zorder = max([_.zorder for _ in ax.get_children()]), cmap = cmap)

    # set yticks
    ax.set_yticks(my_range, df[ylab]);

    # aesthetics
    xlim = ax.get_xlim()
    ax.set_xlim([-np.max(np.abs(xlim)), np.max(np.abs(xlim))])
    if grid_on:
        ax.grid(visible = True, linestyle = '--')
    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
    ax.set_xlabel('NES', fontsize = label_fontsize)
    ax.set_ylabel('Pathways', fontsize = label_fontsize)
    if colorbar:
        cbar = fig.colorbar(im1, shrink = 0.5)
        cbar.set_label(colorbar_label, size = label_fontsize)
        cbar.ax.tick_params(labelsize=tick_fontsize)

    if savefig:
        fig.savefig(savefig, bbox_inches = 'tight', dpi = 150)

    return fig, ax
