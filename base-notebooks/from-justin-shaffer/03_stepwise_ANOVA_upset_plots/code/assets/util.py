import numpy as np
import pandas as pd
from biom import Table
from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance import mantel
from skbio.stats.composition import closure
from numpy.random import poisson, lognormal, normal, randint
from skbio.stats import subsample_counts
from skbio.stats.distance import permanova, mantel
from skbio.stats.distance import DistanceMatrix
from qiime2.plugins.deicode.actions import rpca
from qiime2.plugins.feature_table.actions import rarefy
from qiime2.plugins.diversity.actions import beta_group_significance
from qiime2.plugins.diversity.actions import (beta,
                                              beta_phylogenetic)
from qiime2.plugins.fragment_insertion.actions import filter_features

def mantel_matched(dist, mf_dist, grouping, ids, permutations=5000):
    
    """
    This function creates sample matched
    distance matricies from "grouping"
    and then runs a mantel test.
    
    dist: skbio.stats.distance.DistanceMatrix
    
    mf_dist:
    
    grouping:
    
    ids:
    
    
    """
    
    # check grouping is only == 2
    if len(set(mf_dist[grouping])) != 2:
        raise ValueError('Grouping must have'
                         'exactly two groups.')
    # get org. index name
    ind_name = mf_dist.index.name
    mf_dist = mf_dist.dropna(subset=[ids])
    # get shared IDs
    ids_index = [list(v[ids]) for _, v in mf_dist.groupby(grouping)]
    ids_index = list(set(ids_index[0]).intersection(*ids_index))
    # get groupings and set index as ID
    # & reindex the groupings by same subjects
    grpdf = {k:df.reset_index().set_index(ids).reindex(ids_index)
             for k,df in mf_dist.groupby(grouping)}
    # make distance to compare
    dists = []
    for group_, mf_ in grpdf.items():
        dist_tmp = dist.copy()
        # filter to be matched
        dist_tmp = dist_tmp.filter(mf_[ind_name].values)
        # rename to be matched
        dist_tmp = dist_tmp.to_data_frame()
        rename_ = dict(mf_dist.loc[mf_[ind_name], ids])
        dist_tmp = dist_tmp.rename(rename_,axis=1).rename(rename_,axis=0)
        dist_tmp = DistanceMatrix(dist_tmp.values, dist_tmp.index)
        # save
        dists.append(dist_tmp)
    # run mantel
    return mantel(dists[0], dists[1],
                  permutations=permutations)

def simulate_depth(btsub, depth, read_std):
    
    """
    This function will create a simultion
    of read depth +/- the read std
    on a real table using a Poisson - Log-Normal
    model of counts on a subsampled table
    """
    
    # check depth
    if depth > btsub.sum().min():
        raise ValueError('depth must be '
                         '> min(sum(count))')
    # subsample
    btsub = pd.DataFrame([subsample_counts(btsub.loc[:,s].astype(int), depth)
                          for s in btsub.columns],
                          btsub.columns, btsub.index).T
    btsub = btsub[btsub.sum(axis=1) > 0]
    # Poisson - Log-Normal Counts
    # to get realistic distribution
    sim = closure(btsub.values.T)
    mean = depth
    stdev = read_std * 2
    phi = (stdev ** 2 \
           + mean ** 2) ** 0.5
    mu = mean**2 / phi
    mu = np.log(mu * sim.T)
    sigma = (np.log(phi ** 2 / mean ** 2)) ** 0.5
    sim = np.vstack([poisson(lognormal(mu[:, i],
                                       sigma))
                     for i in range(sim.shape[0])])
    # build table and deouble check
    btsub = pd.DataFrame(sim,
                         btsub.columns,
                         btsub.index).T
    btsub = btsub[btsub.sum(axis=1) > 0]
    
    return btsub

def all_dists(table, rare_depth, tree, minf=None, rpca_depth=None):
    
    """
    Returns basic betas neeeded
    needs both table and rare-depth
    in Q2 artifact format.
    """
    
    table_rare = rarefy(table, rare_depth).rarefied_table
    beta_res = {}
    beta_res['Jaccard'] = beta(table_rare, 'jaccard')
    beta_res['Unweighted UniFrac'] = beta_phylogenetic(table_rare, tree,
                                                       'unweighted_unifrac')    
    beta_res['Weighted UniFrac'] = beta_phylogenetic(table_rare, tree,
                                                       'weighted_unifrac') 
    if rpca_depth == None:
        rpca_depth = rare_depth #use raredepth
    if minf == None:
        minf = 10 #default
    beta_res['RPCA'] = rpca(table, min_feature_count=0,
                            min_sample_count=rare_depth)
    
    return beta_res

def nested_permanova(distdpth, mfdpth, grouping, evaluation, permutations=5000):
    
    """
    This function will run permanova within
    two groups. This is _not_ pairwise!
    """

    # match tables
    shared_ = set(distdpth.ids) & set(mfdpth.index)
    mfdpth = mfdpth.reindex(shared_)
    distdpth = distdpth.filter(shared_)
    # check
    if sorted(mfdpth.index) != sorted(distdpth.ids):
        raise ValueError('Dist and mapping'
                         'could not be matched')
    # group
    perm_res = {}
    for g1_, g1df_ in mfdpth.groupby(grouping):
        # subset Distance
        distdpth_g1 = distdpth.copy()
        distdpth_g1 = distdpth_g1.filter(g1df_.index)
        # run permanova
        perm_res[g1_] = permanova(distdpth_g1,
                                  mfdpth.loc[distdpth_g1.ids,
                                             evaluation],
                                  permutations=permutations)
    return pd.DataFrame(perm_res)

    
def add_taxsplit(taxdf):
    # split taxonomy 
    def tax_split(tax_id, tax_level): return tax_id.split(
        tax_level)[1].split(';')[0]
    for level, lname in zip(['k__', 'p__', 'c__', 'o__',
                             'f__', 'g__', 's__'],
                            ['kingdom', 'phylum', 'class',
                             'order', 'family', 'genus',
                             'species']):
        if lname not in taxdf.columns:
            taxonomy_tmp = []
            for tax in taxdf.Taxon:
                if tax is not np.nan and\
                   level in tax and\
                   len(tax_split(tax, level)) > 0:
                    taxonomy_tmp.append(tax_split(tax,
                                                  level))
                else:
                    taxonomy_tmp.append(np.nan)
            taxdf[lname] = taxonomy_tmp
    return taxdf