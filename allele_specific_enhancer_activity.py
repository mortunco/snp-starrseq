import numpy as np
import pandas as pd
import os
from statsmodels.discrete.discrete_model import NegativeBinomial as NB
from optparse import OptionParser


def merge_rna_dna(var, barcode_allele, rna, dna, allele):
    seq = barcode_allele.start_aln_UMI[(barcode_allele.mutation == var) & (barcode_allele.type == allele)].to_list()
    seq = set(seq) & set(dna.index.to_list()) & set(rna.index.to_list())
    df = pd.DataFrame({'rna': rna.loc[seq].rna_read_count.values, 
                       'dna': dna.loc[seq].dna_read_count.values, 
                       'type':0}, 
                        index = seq)
    return df

def parse_variant(s):
    return int(s.split(';')[0].replace('chr','')), int(s.split(';')[1])

def merge_rna_dna(var, barcode_allele, rna, dna, allele):
    seq = barcode_allele.start_aln_UMI[(barcode_allele.mutation == var) & (barcode_allele.type == allele)].to_list()
    seq = set(seq) & set(dna.index.to_list()) & set(rna.index.to_list())
    df = pd.DataFrame({'rna': rna.loc[seq].rna_read_count.values, 
                       'dna': dna.loc[seq].dna_read_count.values, 
                       'type':allele, 'intercept':1}, 
                        index = seq)
    return df
def NBregression(analysis_df):
    try:
        fit = NB(endog = analysis_df['rna'].values, 
                             exog = analysis_df[['type', 'intercept']].values,
                             loglike_method='nb2', 
                             offset=np.log(analysis_df['dna'].values)).fit()

        pvalue = fit.pvalues[0]
        coef = fit.params[0]
        return coef, pvalue
    except:
        return np.nan, np.nan

def NB_for_one_variant(var, rna, dna, barcode_allele):
    chrom, pos = parse_variant(var)
    ref_df = merge_rna_dna(var, barcode_allele, rna, dna,1)
    alt_df = merge_rna_dna(var, barcode_allele, rna, dna,2)
    ref_count = len(ref_df)
    alt_count = len(alt_df)

    if (ref_count >= 3) & (alt_count >= 3):
        analysis_df = pd.concat([ref_df, alt_df], axis = 0)
        coef, pvalue =  NBregression(analysis_df)
    else:
        coef = pvalue = np.nan
    output = pd.Series({'var': var, 'chr': chrom, 'pos': pos, 'log2fc': np.log2(np.exp(coef)), 'pvalue': pvalue, 
           'ref_count': len(ref_df), 'alt_count': len(alt_df)})
    return output

def NB_for_all_variant(var_list, rna, dna, barcode_allele):
    return pd.DataFrame(NB_for_one_variant(var, rna, dna, barcode_allele) for var in var_list)
    
parser = OptionParser()
parser.add_option("--dna", dest="dna", metavar="FILE")
parser.add_option("--rna", dest="rna", metavar="FILE")
parser.add_option("--barcode_allele", dest="barcode_allele", metavar="FILE")
(options, args) = parser.parse_args()



rna = pd.read_csv(options.rna, sep = ' ', header = None, names = ['rna_read_count', 'start_aln_UMI', 'fragment_name'] )
rna = rna[['start_aln_UMI', 'rna_read_count']].groupby('start_aln_UMI').sum()
dna = pd.read_csv(options.dna, sep = ' ', header = None, names = ['dna_read_count', 'start_aln_UMI', 'fragment_name'])
dna = dna[['start_aln_UMI', 'dna_read_count']].groupby('start_aln_UMI').sum()
barcode_allele = pd.read_csv(options.barcode_allele, sep = '\t')

var_list = list(barcode_allele.mutation.unique())


output = NB_for_all_variant(var_list, rna, dna, barcode_allele)
output.to_csv(os.path.splitext(options.rna)[0] + '.nb', index = False, sep = '\t')

