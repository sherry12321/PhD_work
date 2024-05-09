import numpy as np
import scanpy as sc
import pandas as pd
import os
from typing import Optional, Union
import itertools
from anndata import AnnData
from celltypist import logger
from celltypist.models import Model
from scipy.sparse import spmatrix
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import StandardScaler
from sklearn import __version__ as skv
from datetime import datetime
import sys
sys.path.append(os.path.abspath(os.getcwd()))
import classifier


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

def _SGDClassifier(indata, labels,
                   alpha, max_iter, n_jobs,
                   mini_batch, batch_number, batch_size, epochs, balance_cell_type, penalty ,  **kwargs) -> SGDClassifier:
    """
    For internal use 
    
    ONE NEW ARG
    penalty
        allows to user specify what type of regularization
    """
    loss_mode = 'log_loss' if float(skv[:3]) >= 1.1 else 'log'
    classifier = SGDClassifier(loss = loss_mode, penalty = penalty, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, **kwargs)
    if not mini_batch:
        logger.info(f"🏋️ Training data using SGD logistic regression")
        #print('')
        if (len(labels) > 100000) and (indata.shape[1] > 10000):            
            logger.warn(f"⚠️ Warning: it may take a long time to train this dataset with {len(labels)} cells and {indata.shape[1]} genes, try to downsample cells and/or restrict genes to a subset (e.g., hvgs)")
        classifier.fit(indata, labels)
    else:
        logger.info(f"🏋️ Training data using mini-batch SGD logistic regression")
        #print('')
        no_cells = len(labels)
        if no_cells < 10000:
            logger.warn(f"⚠️ Warning: the number of cells ({no_cells}) is not big enough to conduct a proper mini-batch training. You may consider using traditional SGD classifier (mini_batch = False)")
            #print('')
        if no_cells <= batch_size:
            raise ValueError(
                    f"🛑 Number of cells ({no_cells}) is fewer than the batch size ({batch_size}). Decrease `batch_size`, or use SGD directly (mini_batch = False)")
        no_cells_sample = min([batch_number*batch_size, no_cells])
        starts = np.arange(0, no_cells_sample, batch_size)
        if balance_cell_type:
            celltype_freq = np.unique(labels, return_counts = True)
            len_celltype = len(celltype_freq[0])
            mapping = pd.Series(1 / (celltype_freq[1]*len_celltype), index = celltype_freq[0])
            p = mapping[labels].values
        for epoch in range(1, (epochs+1)):
            logger.info(f"⏳ Epochs: [{epoch}/{epochs}]")
            #print('')
            if not balance_cell_type:
                sampled_cell_index = np.random.choice(no_cells, no_cells_sample, replace = False)
            else:
                sampled_cell_index = np.random.choice(no_cells, no_cells_sample, replace = False, p = p)
            for start in starts:
                classifier.partial_fit(indata[sampled_cell_index[start:start+batch_size]], labels[sampled_cell_index[start:start+batch_size]], classes = np.unique(labels))
    return classifier

def  celltypist_train_modified(X = None,
          labels: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
          genes: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
          transpose_input: bool = False,
          with_mean: bool = True,
          check_expression: bool = False,
          normalize: bool = False,
          #LR param
          C: float = 1.0, solver: Optional[str] = None, max_iter: Optional[int] = None, n_jobs: Optional[int] = None,
          #SGD param
          use_SGD: bool = False, alpha: float = 0.0001,
          #mini-batch
          mini_batch: bool = False, batch_number: int = 100, batch_size: int = 1000, epochs: int = 10, balance_cell_type: bool = False,
          #feature selection
          feature_selection: bool = False, top_genes: int = 300, use_additional: bool = False, additional_genes: Optional[np.ndarray] = None,
          #description
          date: str = '', details: str = '', url: str = '', source: str = '', version: str = '',
          #penalty
          penalty: str = 'l2', switch_penalty: bool = False,
          #other param
          **kwargs): 
    """
    Train a celltypist model using mini-batch (optional) logistic classifier with a global solver or stochastic gradient descent (SGD) learning. 
    A version of the celltypist fxn train that adds some additional choices/functions 

    Parameters
    ----------
    X
        Path to the input count matrix (supported types are csv, txt, tsv, tab and mtx) or AnnData (h5ad).
        Also accepts the input as an :class:`~anndata.AnnData` object, or any array-like objects already loaded in memory.
        See `check_expression` for detailed format requirements.
        A cell-by-gene format is desirable (see `transpose_input` for more information).
    labels
        Path to the file containing cell type label per line corresponding to the cells in `X`.
        Also accepts any list-like objects already loaded in memory (such as an array).
        If `X` is specified as an AnnData, this argument can also be set as a column name from cell metadata.
    genes
        Path to the file containing one gene per line corresponding to the genes in `X`.
        Also accepts any list-like objects already loaded in memory (such as an array).
        Note `genes` will be extracted from `X` where possible (e.g., `X` is an AnnData or data frame).
    transpose_input
        Whether to transpose the input matrix. Set to `True` if `X` is provided in a gene-by-cell format.
        (Default: `False`)
    with_mean
        Whether to subtract the mean values during data scaling. Setting to `False` can lower the memory usage when the input is a sparse matrix but may slightly reduce the model performance.
        (Default: `True`)
    check_expression
        Check whether the expression matrix in the input data is supplied as required.
        Except the case where a path to the raw count table file is specified, all other inputs for `X` should be in log1p normalized expression to 10000 counts per cell.
        Set to `False` if you want to train the data regardless of the expression formats.
        (Default: `True`)
    C
        Inverse of L2 regularization strength for traditional logistic classifier. A smaller value can possibly improve model generalization while at the cost of decreased accuracy.
        This argument is ignored if SGD learning is enabled (`use_SGD = True`).
        (Default: 1.0)
    solver
        Algorithm to use in the optimization problem for traditional logistic classifier.
        The default behavior is to choose the solver according to the size of the input data.
        This argument is ignored if SGD learning is enabled (`use_SGD = True`).
    max_iter
        Maximum number of iterations before reaching the minimum of the cost function.
        Try to decrease `max_iter` if the cost function does not converge for a long time.
        This argument is for both traditional and SGD logistic classifiers, and will be ignored if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        Default to 200, 500, and 1000 for large (>500k cells), medium (50-500k), and small (<50k) datasets, respectively.
    n_jobs
        Number of CPUs used. Default to one CPU. `-1` means all CPUs are used.
        This argument is for both traditional and SGD logistic classifiers.
    use_SGD
        Whether to implement SGD learning for the logistic classifier.
        (Default: `False`)
    alpha
        L2 regularization strength for SGD logistic classifier. A larger value can possibly improve model generalization while at the cost of decreased accuracy.
        This argument is ignored if SGD learning is disabled (`use_SGD = False`).
        (Default: 0.0001)
    mini_batch
        Whether to implement mini-batch training for the SGD logistic classifier.
        Setting to `True` may improve the training efficiency for large datasets (for example, >100k cells).
        This argument is ignored if SGD learning is disabled (`use_SGD = False`).
        (Default: `False`)
    batch_number
        The number of batches used for training in each epoch. Each batch contains `batch_size` cells.
        For datasets which cannot be binned into `batch_number` batches, all batches will be used.
        This argument is relevant only if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        (Default: 100)
    batch_size
        The number of cells within each batch.
        This argument is relevant only if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        (Default: 1000)
    epochs
        The number of epochs for the mini-batch training procedure.
        The default values of `batch_number`, `batch_size`, and `epochs` together allow observing ~10^6 training cells.
        This argument is relevant only if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        (Default: 10)
    balance_cell_type
        Whether to balance the cell type frequencies in mini-batches during each epoch.
        Setting to `True` will sample rare cell types with a higher probability, ensuring close-to-even cell type distributions in mini-batches.
        This argument is relevant only if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        (Default: `False`)
    feature_selection
        Whether to perform two-pass data training where the first round is used for selecting important features/genes using SGD learning.
        If `True`, the training time will be longer.
        (Default: `False`)
    top_genes
        The number of top genes selected from each class/cell-type based on their absolute regression coefficients.
        The final feature set is combined across all classes (i.e., union).
        (Default: 300)
    date
        Free text of the date of the model. Default to the time when the training is completed.
    details
        Free text of the description of the model.
    url
        Free text of the (possible) download url of the model.
    source
        Free text of the source (publication, database, etc.) of the model.
    version
        Free text of the version of the model.
    **kwargs
        Other keyword arguments passed to :class:`~sklearn.linear_model.LogisticRegression` (`use_SGD = False`) or :class:`~sklearn.linear_model.SGDClassifier` (`use_SGD = True`).
        
    FIVE NEW PARAMETERS: 
    penalty
        Which regularization method to use 
        (Default: "l2")
    switch_penalty
        Whether to switch the type of regualarization for the second round of model training 
        (Default: False)
    use_additional 
        Whether to confirm that additional genes are included in feature_selection (they are added if not)
        This argument is relevant only if feature selection happens (`feature_selection = True`) 
        (Default: False)
    additional_genes
        List of gene names from ctyopus cell identities dictionary
        This argument is relevant only if feature selection with additional genes happens (`feature_selection = True` and `use_additional = True`)
    normalize 
        Whether or not to normalize the data in prepare_data(). If True, data will be normalized to median library size. Will be ignored if X is not given as AnnData. 
        (Default: `False`)

    Returns
    ----------
    :class:`~celltypist.models.Model`
        An instance of the :class:`~celltypist.models.Model` trained by celltypist.
    pandas.DataFrame
        A DataFrame that contains the 300 genes that were selected for each cell type (organized by cell type)
    """
    #prepare
    logger.warn("🍳 Preparing data before training")
    #print('')
    indata, labels, genes = _prepare_data(X, labels, genes, transpose_input, normalize)
    if isinstance(indata, pd.DataFrame):
        indata = indata.values
    elif with_mean and isinstance(indata, spmatrix):
        indata = indata.toarray()
    labels = np.array(labels)
    genes = np.array(genes)
    #check
    ##NEED TO CHANGE 10000 TO MEDIAN AMOUNT 
    print("Check expression: ", check_expression)
    print("sum: ", (np.abs(np.expm1(indata[0]).sum()-10000) > 1))
    if check_expression and (np.abs(np.expm1(indata[0]).sum()-10000) > 1):
        raise ValueError(
                "🛑 Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell")
    if len(labels) != indata.shape[0]:
        raise ValueError(
                f"🛑 Length of training labels ({len(labels)}) does not match the number of input cells ({indata.shape[0]})")
    if len(genes) != indata.shape[1]:
        raise ValueError(
                f"🛑 The number of genes ({len(genes)}) provided does not match the number of genes in the training data ({indata.shape[1]})")
    #filter
    flag = indata.sum(axis = 0) == 0
    if isinstance(flag, np.matrix):
        flag = flag.A1
    if flag.sum() > 0:
        logger.info(f"✂️ {flag.sum()} non-expressed genes are filtered out")
        #print('')
        #indata = indata[:, ~flag]
        genes = genes[~flag]
    #report data stats
    logger.info(f"🔬 Input data has {indata.shape[0]} cells and {(~flag).sum()} genes")
    #print('')
    #scaler
    logger.info(f"⚖️ Scaling input data")
    #print('')
    
    scaler = StandardScaler(with_mean = with_mean)
    indata = scaler.fit_transform(indata[:, ~flag] if flag.sum() > 0 else indata)
    indata[indata > 10] = 10
    #sklearn (Cython) does not support very large sparse matrices for the time being
    if isinstance(indata, spmatrix) and ((indata.indices.dtype == 'int64') or (indata.indptr.dtype == 'int64')):
        indata = indata.toarray()
    #max_iter
    if max_iter is None:
        if indata.shape[0] < 50000:
            max_iter = 1000
        elif indata.shape[0] < 500000:
            max_iter = 500
        else:
            max_iter = 200
    #classifier
    if use_SGD or feature_selection:
        classifier = _SGDClassifier(indata = indata, labels = labels, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, mini_batch = mini_batch, batch_number = batch_number, batch_size = batch_size, epochs = epochs, balance_cell_type = balance_cell_type, penalty = penalty, **kwargs)
    else:
        classifier = _LRClassifier(indata = indata, labels = labels, C = C, solver = solver, max_iter = max_iter, n_jobs = n_jobs, **kwargs)
    #feature selection -> new classifier and scaler
    genes_dict = {}
    if feature_selection:
        logger.info(f"🔎 Selecting features")
        #print('')
        if len(genes) <= top_genes:
            raise ValueError(
                    f"🛑 The number of genes ({len(genes)}) is fewer than the `top_genes` ({top_genes}). Unable to perform feature selection")
        gene_index = np.argpartition(np.abs(classifier.coef_), -top_genes, axis = 1)[:, -top_genes:]
        
        #making a dictionary of all of the genes that are selected for one cell
        #return gene_index
        for row in range(len(gene_index)):
            gene_names = genes[gene_index[row]]
            cell_name = classifier.classes_[row]
            genes_dict[cell_name] = gene_names    
        gene_index = np.unique(gene_index)
        
        if use_additional: 
            logger.info(f"🧬 {len(gene_index)} features are selected before checking with external genes")
            #print('')
            #confirming that all external genes are in the top genes used in feature selection
            #first get a list of all the indexs of additional_genes 
            external_gene_index = []
            for x in additional_genes:
                if x in genes: 
                    idx = np.where(genes==x)[0][0]
                    external_gene_index.append(idx)
            for x in external_gene_index: 
                if x not in gene_index: 
                    gene_index = np.append(gene_index, x)
            gene_index = np.unique(gene_index)
            logger.info(f"🧬 {len(gene_index)} features are selected after checking with external genes")
            #print('')
        else:
            logger.info(f"🧬 {len(gene_index)} features are selected")
            #print('')
        genes = genes[gene_index]
        
        logger.info(f"🏋️ Starting the second round of training")
        #print('')
        if switch_penalty: 
            if penalty == "l2":
                penalty = "l1"
            else: 
                penalty = "l2"
        if use_SGD:
            classifier = _SGDClassifier(indata = indata[:, gene_index], labels = labels, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, mini_batch = mini_batch, batch_number = batch_number, batch_size = batch_size, epochs = epochs, balance_cell_type = balance_cell_type, penalty = penalty, **kwargs)
        else:
            classifier = _LRClassifier(indata = indata[:, gene_index], labels = labels, C = C, solver = solver, max_iter = max_iter, n_jobs = n_jobs, **kwargs)
        scaler.mean_ = scaler.mean_[gene_index]
        scaler.var_ = scaler.var_[gene_index]
        scaler.scale_ = scaler.scale_[gene_index]
        scaler.n_features_in_ = len(gene_index)
    #model finalization
    classifier.features = genes
    classifier.n_features_in_ = len(genes)
    if not date:
        date = str(datetime.now())
    description = {'date': date, 'details': details, 'url': url, 'source': source, 'version': version, 'number_celltypes': len(classifier.classes_)}
    logger.info(f"✅ Model training done!")
    #print('')
    genes_df = pd.DataFrame(genes_dict)
    genes_df = genes_df.T
    return Model(classifier, scaler, description), genes_df


def _prepare_data(X, labels, genes, transpose, normalize) -> tuple:
    """
    For internal use. Prepare data for celltypist training.
    
    ONE NEW ARG
    normalize 
        whether or not to normalize anndata input
    """
    if (X is None) or (labels is None):
        raise Exception(
                "🛑 Missing training data and/or training labels. Please provide both arguments")
    if isinstance(X, AnnData) or (isinstance(X, str) and X.endswith('.h5ad')):
        adata = sc.read(X) if isinstance(X, str) else X
        adata.var_names_make_unique()
        
        if normalize:
            adata.X = adata.raw.X
            adata.X= np.log1p(adata.X)
            sc.pp.normalize_total(adata)
            indata = adata.X
        else: 
            indata = adata.X
            
        if adata.X.min() < 0:
            logger.info("👀 Detected scaled expression in the input data, will try the .raw attribute")
            #print('')
            try:
                indata = adata.raw.X
                genes = adata.raw.var_names
            except Exception as e:
                raise Exception(
                        f"🛑 Fail to use the .raw attribute in the input object. {e}")
        else:
            #indata = adata.X
            genes = adata.var_names
        if isinstance(labels, str) and (labels in adata.obs):
            labels = adata.obs[labels]
        else:
            labels = _to_vector(labels)
    elif isinstance(X, str) and X.endswith(('.csv', '.txt', '.tsv', '.tab', '.mtx', '.mtx.gz')):
        adata = sc.read(X)
        if transpose:
            adata = adata.transpose()
        if X.endswith(('.mtx', '.mtx.gz')):
            if genes is None:
                raise Exception(
                        "🛑 Missing `genes`. Please provide this argument together with the input mtx file")
            genes = _to_vector(genes)
            if len(genes) != adata.n_vars:
                raise ValueError(
                        f"🛑 The number of genes provided does not match the number of genes in {X}")
            adata.var_names = np.array(genes)
        adata.var_names_make_unique()
        if not float(adata.X.max()).is_integer():
            logger.warn(f"⚠️ Warning: the input file seems not a raw count matrix. The trained model may be biased")
            #print('')
        sc.pp.normalize_total(adata, target_sum=1e4) 
        sc.pp.log1p(adata)
        indata = adata.X
        genes = adata.var_names
        labels = _to_vector(labels)
    elif isinstance(X, str):
        raise ValueError(
                "🛑 Invalid input. Supported types: .csv, .txt, .tsv, .tab, .mtx, .mtx.gz and .h5ad")
    else:
        logger.info("👀 The input training data is processed as an array-like object")
        #print('')
        indata = X
        if transpose:
            indata = indata.transpose()
        if isinstance(indata, pd.DataFrame):
            genes = indata.columns
        else:
            if genes is None:
                raise Exception(
                        "🛑 Missing `genes`. Please provide this argument together with the input training data")
            genes = _to_vector(genes)
        labels = _to_vector(labels)
    return indata, labels, genes

def _to_vector(_vector_or_file):
    """
    For internal use. Turn a file into an array.
    """
    if isinstance(_vector_or_file, str):
        try:
            return pd.read_csv(_vector_or_file, header=None)[0].values
        except Exception as e:
            raise Exception(
                    f"🛑 {e}")
    else:
        return _vector_or_file




def train_test_split(adata, frac: int = 0.75):
    """
    USING OUTLINE OF CODE FROM trVAE https://doi.org/10.1093/bioinformatics/btaa800
    Split AnnData into test and train datasets - maintains annotations
    
    Params: 
    adata
        Annotated data matrix (Anndata)
    frac
        Fraction of cells to be used in the training set
    """
    no_idx_train = int(adata.shape[0] * frac)
    indices = np.arange(adata.shape[0])
    np.random.shuffle(indices)
    train_index = indices[:no_idx_train]
    test_index = indices[no_idx_train:]
    train = adata[train_index]
    test = adata[test_index]
    return train, test



def annotate(filename: Union[AnnData,str] = "",
             model: Optional[Union[str, Model]] = None,
             transpose_input: bool = False,
             gene_file: Optional[str] = None,
             cell_file: Optional[str] = None,
             mode: str = 'best match',
             p_thres: float = 0.5,
             majority_voting: bool = False,
             over_clustering: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
             min_prop: float = 0, 
             check_expression: bool = False) -> classifier.AnnotationResult:
    """
    Run the prediction and (optional) majority voting to annotate the input dataset.

    Parameters
    ----------
    filename
        Path to the input count matrix (supported types are csv, txt, tsv, tab and mtx) or AnnData (h5ad).
        If it's the former, a cell-by-gene format is desirable (see `transpose_input` for more information).
        Also accepts the input as an :class:`~anndata.AnnData` object already loaded in memory.
        Genes should be gene symbols. Non-expressed genes are preferred to be provided as well.
    model
        Model used to predict the input cells. Default to using the `'Immune_All_Low.pkl'` model.
        Can be a :class:`~celltypist.models.Model` object that wraps the logistic Classifier and the StandardScaler, the
        path to the desired model file, or the model name.
        To see all available models and their descriptions, use :func:`~celltypist.models.models_description`.
    transpose_input
        Whether to transpose the input matrix. Set to `True` if `filename` is provided in a gene-by-cell format.
        (Default: `False`)
    gene_file
        Path to the file which stores each gene per line corresponding to the genes used in the provided mtx file.
        Ignored if `filename` is not provided in the mtx format.
    cell_file
        Path to the file which stores each cell per line corresponding to the cells used in the provided mtx file.
        Ignored if `filename` is not provided in the mtx format.
    mode
        The way cell prediction is performed.
        For each query cell, the default (`'best match'`) is to choose the cell type with the largest score/probability as the final prediction.
        Setting to `'prob match'` will enable a multi-label classification, which assigns 0 (i.e., unassigned), 1, or >=2 cell type labels to each query cell.
        (Default: `'best match'`)
    p_thres
        Probability threshold for the multi-label classification. Ignored if `mode` is `'best match'`.
        (Default: 0.5)
    majority_voting
        Whether to refine the predicted labels by running the majority voting classifier after over-clustering.
        (Default: `False`)
    over_clustering
        This argument can be provided in several ways:
        1) an input plain file with the over-clustering result of one cell per line.
        2) a string key specifying an existing metadata column in the AnnData (pre-created by the user).
        3) a python list, tuple, numpy array, pandas series or index representing the over-clustering result of the input cells.
        4) if none of the above is provided, will use a heuristic over-clustering approach according to the size of input data.
        Ignored if `majority_voting` is set to `False`.
    min_prop
        For the dominant cell type within a subcluster, the minimum proportion of cells required to support naming of the subcluster by this cell type.
        Ignored if `majority_voting` is set to `False`.
        Subcluster that fails to pass this proportion threshold will be assigned `'Heterogeneous'`.
        (Default: 0)

    Returns
    ----------
    :class:`~celltypist.classifier.AnnotationResult`
        An :class:`~celltypist.classifier.AnnotationResult` object. Four important attributes within this class are:
        1) :attr:`~celltypist.classifier.AnnotationResult.predicted_labels`, predicted labels from celltypist.
        2) :attr:`~celltypist.classifier.AnnotationResult.decision_matrix`, decision matrix from celltypist.
        3) :attr:`~celltypist.classifier.AnnotationResult.probability_matrix`, probability matrix from celltypist.
        4) :attr:`~celltypist.classifier.AnnotationResult.adata`, AnnData representation of the input data.
    """
    #load model
    lr_classifier = model if isinstance(model, Model) else Model.load(model)
    #construct Classifier class
    clf = classifier.Classifier(filename = filename, model = lr_classifier, transpose = transpose_input, gene_file = gene_file, cell_file = cell_file)
    #predict
    predictions = clf.celltype(mode = mode, p_thres = p_thres)
    if not majority_voting:
        return predictions
    if predictions.cell_count <= 50:
        logger.warn(f"⚠️ Warning: the input number of cells ({predictions.cell_count}) is too few to conduct proper over-clustering; no majority voting is performed")
        return predictions
    #over clustering
    if over_clustering is None:
        over_clustering = clf.over_cluster()
        predictions.adata = clf.adata
    elif isinstance(over_clustering, str):
        if over_clustering in clf.adata.obs:
            over_clustering = clf.adata.obs[over_clustering]
        else:
            logger.info(f"👀 Did not identify '{over_clustering}' as a cell metadata column, assume it to be a plain text file")
            try:
                with open(over_clustering, 'rt') as f:
                    over_clustering = [x.strip() for x in f.readlines()]
            except Exception as e:
                raise Exception(
                        f"🛑 {e}")
    if len(over_clustering) != clf.adata.n_obs:
        raise ValueError(
                f"🛑 Length of `over_clustering` ({len(over_clustering)}) does not match the number of input cells ({clf.adata.n_obs})")
    #majority voting
    return classifier.Classifier.majority_vote(predictions, over_clustering, min_prop = min_prop)

