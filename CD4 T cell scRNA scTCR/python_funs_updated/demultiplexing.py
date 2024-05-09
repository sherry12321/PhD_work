import joypy
import anndata
import demuxEM
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp_sparse
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from skimage import filters
from scipy.signal import find_peaks
from scipy.optimize import *


def asinh_trans(_ann):
    tmpm=np.random.normal(loc=0.5, scale=0.5, size=_ann.X.shape)
    x_trans = tmpm + _ann.X.toarray()
    _ann.X = np.arcsinh(x_trans - min(x_trans.flatten()))
    
    
def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2))
    
    
def find_thre(ab_exp, hash_num):
    # get fit and find peak
    bw=0.1
    kde = KernelDensity(bandwidth=bw, kernel='gaussian')
    kde.fit(ab_exp.flatten().reshape(-1,1))
    _x = np.linspace(np.min(ab_exp), np.max(ab_exp), 200).reshape(-1,1)
    x_fit = np.exp(kde.score_samples(_x))
    x_fit -= min(x_fit)
    x_fit /= max(x_fit)
    peak,ph = find_peaks(x_fit,height=1./hash_num/5.,width=1.5,distance=20)
    
    if len(peak)<2:
        # unimodal distribution, use Gaussian fit
        n_bins = min(50,int(len(ab_exp)/100))
        counts, edges = np.histogram(ab_exp,n_bins,density=True)
        bin_center = edges[:-1] + (edges[1] - edges[0])/2.
        popt_gauss, pcov_gauss = curve_fit(gaussian, bin_center, counts, p0=[0.6,_x[peak[0]][0],2])
        thre_lower = popt_gauss[1] - 3.09*abs(popt_gauss[2])
        thre_upper = popt_gauss[1] + 3.09*abs(popt_gauss[2])
        return max(thre_upper,thre_lower)
    
    # multimodal distribution, find the lowest valley
    thre_final = 0
    valley_counts_ratio = 99
    for i in range(len(peak)-1):
        vally_bound_low = _x[peak[i]]
        vally_bound_high = _x[peak[i+1]]
        vally_region = ab_exp[np.logical_and(ab_exp>vally_bound_low,ab_exp<vally_bound_high)]
        try:
            thre = filters.threshold_minimum(vally_region)
        except:
            thre = filters.threshold_minimum(ab_exp)
    
        thre_arg = np.argmin((_x.reshape(-1,1)-thre)**2)
        ratio_temp = x_fit[thre_arg]/np.min([x_fit[peak[t]] for t in range(len(peak))])
        if ratio_temp < valley_counts_ratio:
            thre_final = thre
            valley_counts_ratio = ratio_temp
    return thre_final


def get_demux(_ann,method='demuxEM'):
    # generate RNA and HTO matrix
    ann_rna = _ann[:, _ann.var['feature_types']=='Gene Expression'].copy()
    if not ann_rna.n_obs:
        ann_rna = _ann[:,['ADT' in t for t in _ann.var['gene_ids']]].copy()
    ann_hto = _ann[:, ['HTO' in t for t in _ann.var['gene_ids']]].copy()

    # filter out 0 HTO counts cells for cleaner demuxEM input
    umi_hto = np.squeeze(np.asarray(ann_hto.X.sum(axis=1)))
    demux_filt = umi_hto > 0
    ann_rna = ann_rna[demux_filt].copy()
    ann_hto = ann_hto[demux_filt].copy()
    
    if method == 'demuxEM':
        # use demuxEM to get demultiplexing results
        demuxEM.estimate_background_probs(ann_hto)
        demuxEM.demultiplex(ann_rna, ann_hto)

        # clean up assignment, replace detailed assignment list to doublets and negatives
        ann_rna.obs['assignment'] = ann_rna.obs['assignment'].tolist()
        ann_rna.obs.loc[ann_rna.obs['demux_type'] == 'doublet', 'assignment'] = 'Doublet'
        ann_rna.obs.loc[ann_rna.obs['demux_type'] == 'unknown', 'assignment'] = 'Negative'

        _ann.obs['assignment'] = ['Negative'] * _ann.n_obs
        if np.sum(umi_hto==0) > 0:
            _ann.obs.loc[ann_rna.obs_names, 'assignment'] = ann_rna.obs['assignment']
        else:
            _ann.obs['assignment'] = ann_rna.obs.loc[_ann.obs_names, 'assignment']
            
    else:
        if 'thre' not in _ann.var.columns:
            _ann.var['thre'] = np.zeros(_ann.n_vars)
        asinh_trans(ann_hto)
        n_hash = ann_hto.n_vars
        ann_hto.obs['assignment'] = ['Negative',]*ann_hto.n_obs
        ann_hto.obs['positive_counts'] = np.zeros(ann_hto.n_obs)
        fig_row = int((ann_hto.n_vars+3)/4)
        fig,ax = plt.subplots(fig_row,4,figsize=(12,fig_row*3))
        ax = ax.flatten()
        ii=0
        for tag in ann_hto.var_names:
            ab = ann_hto.obs_vector(tag)
            ax[ii].hist(ab,100)
            ax[ii].set_title(tag)
            if method == 'update':
                hash_thre = _ann.var.loc[tag, 'thre']
            else:
                hash_thre = find_thre(ab,n_hash)
                _ann.var.loc[tag, 'thre'] = hash_thre
            ax[ii].axvline(hash_thre,c='r')
            ab_posi = ab > hash_thre
            ann_hto.obs.loc[ab_posi, 'assignment'] = tag
            ann_hto.obs.loc[ab_posi, 'positive_counts'] += 1
            ii += 1
        for hist_ax in ax[ii:]:
            hist_ax.axis('off')
        plt.tight_layout()
        ann_hto.obs.loc[ann_hto.obs['positive_counts'] > 1, 'assignment'] = 'Doublets'
        _ann.obs['assignment'] = ['Negative',]*_ann.n_obs
        _ann.obs.loc[ann_hto.obs_names, 'assignment'] = ann_hto.obs['assignment']

        
def get_ridge(_ann):
    df = pd.DataFrame(_ann.X,columns=_ann.var_names.tolist())
    df['demuxem_hashID'] = _ann.obs['assignment'].tolist()
    df = df.sort_values('demuxem_hashID')
    for tags in _ann.var_names:
        fig, axes = joypy.joyplot(df,by='demuxem_hashID',ylim='own',column=tags,title=tags,
                                alpha=0.8,overlap=1,colormap= matplotlib.cm.tab10,figsize=(6,6))
        #plt.savefig(os.path.join(fig_path, 'ridge_'+tags+'_ADT.png'),dpi=200)