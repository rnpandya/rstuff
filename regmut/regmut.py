#
# regmut.py
#
# regulatory mutation analysis support code
#
# (c) 2015 Ravi Pandya, Microsoft Research
#

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import re
import itertools
import os
import os.path
from Bio import motifs
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

# DataFrame of ENCODE TF binding sites
def read_encode_tfbs(dir='d:/sequence'):
    return pd.read_table(os.path.join(dir, 'encode/wgEncodeRegTfbsClusteredV3.bed.gz'),
        header=None,names=['seqname','start','end','gene','score', 'nexpt', 'expts', 'expscores'],
        dtype={'start':np.int, 'end':np.int})

def splitrow(row, c1, c2, sep):
    rows = []
    for s1, s2 in zip(row[c1].split(sep), row[c2].split(sep)):
        row2 = row.to_dict()
        row2[c1] = s1
        row2[c2] = s2
        rows.append(row2)
    return rows

# dict of cancer -> TFBS DataFrame
def read_tfbs_bycancer(tfbs, dir='d:/sequence'):
    x = tfbs
    x = pd.DataFrame(list(itertools.chain(*x.apply(lambda r: splitrow(r, 'expts', 'expscores', ','), axis=1))))
    x = x.rename(columns={'expts': 'expt', 'expscores': 'expscore'})
    x['expt'] = x['expt'].astype(np.int32)
    tfbsx = x
    # read multi-stage mapping from experiment->cell->tissue->tissue type->cancer
    tfmeta = pd.read_csv(os.path.join(dir, 'encode/wgEncodeRegTfbsClusteredInputsV3.tab.gz'),
        sep='\t', names=['id','info','tf','antibody','cell','unk','site'])
    tfmeta = tfmeta.reset_index().rename(columns={'index':'expt'})
    celltypes = pd.read_csv('data/EncodeCellTypes.csv')
    # for editing ... celltypes.Tissue.drop_duplicates().to_csv('EncodeTissues.csv',index=False,header=False)
    celltissue = pd.read_csv('data/EncodeTissues.csv')
    tcgatypes = pd.melt(pd.read_csv('data/TcgaTissues.csv'), id_vars=('TissueType'), var_name='cancer')
    tcgatypes = tcgatypes[tcgatypes.value == 1][['TissueType','cancer']]
    # join them all!
    x = tfbsx.merge(tfmeta[['expt','cell']], how='left', on='expt')
    x = x.rename(columns={'cell':'Cell'}).merge(celltypes[['Cell','Tissue']], how='left', on='Cell')
    x = x.merge(celltissue, how='left', on='Tissue')
    x = x.merge(tcgatypes, how='left', on='TissueType')
    x = {g: tfs[['seqname','start','end','gene']].drop_duplicates() for g, tfs in x.groupby('cancer')}
    x['CRC'] = x['COAD'] # synonym
    return x

# DataFrame of 505 x 2 matched donor/random mutation sets
def read_tcga505_dna_all(dir='d:/sequence'):
    dna = pd.read_table(os.path.join(dir, 'pancan/tcga505/mutations.tsv.gz'), dtype={'chr': np.character})
    dna_rnd = pd.read_table(os.path.join(dir, 'pancan/tcga505/randomised_TCGA505_toshare.txt.gz'),
        header=None,names=['barcode', 'cancer', 'chr', 'pos', 'ref_allele', 'var_allele'],
        dtype={'chr': np.character}, skiprows=1)
    dna['rnd'] = False
    dna_rnd['rnd'] = True
    dna_all = dna.append(dna_rnd)
    dna_all['chr'] = ['chr' + x for x in dna_all['chr']]
    return dna_all

# TF gene symbol -> normalized log odds PWM
def read_jaspar_pwms(file='jaspar/pfm_vertebrates.txt', dir='d:/sequence'):
    with open(os.path.join(dir, file), 'r') as h:
        jaspar = {m.name.strip().split(':')[0].upper(): m.counts.normalize().log_odds() for m in motifs.parse(h, 'jaspar')}
        return jaspar

# reduce tfbs to those for which we have matrices
def filter_tfbs_in(tfbs, tfbs_bycancer, pwms):
    common = set(tfbs['gene']).intersection(set(pwms.keys()))
    tfbs2 = tfbs[tfbs.gene.isin(common)]
    tfbs_bycancer2 = {c: tf[tf.gene.isin(common)] for c,tf in tfbs_bycancer.iteritems()} if tfbs_bycancer != None else None
    return tfbs2, tfbs_bycancer2

# get reference genome as chrname -> Seq
def read_reference(path = 'genomes/hg19.fa', dir='d:/sequence'):
    return {contig.name: Seq(str(contig.seq), IUPAC.unambiguous_dna)
        for contig in SeqIO.parse(os.path.join(dir, path), 'fasta')}

# left join by genomic regions
# cols[12] are names of chrom, start, end followed by names of columns to include in join result
# limitations:
# genomic range columns (rcols1, rcols2) must be distinct
# only works if region sets are non-overlapping
# drops duplicates
def merge_regions(df1, cols1, df2, cols2):
    result = None
    [seqname1, start1, end1] = cols1[:3]
    [seqname2, start2, end2] = cols2[:3]
    cols1 = list(pd.Series(cols1).drop_duplicates())
    cols2 = list(pd.Series(cols2).drop_duplicates())
    for contig in set(df1[seqname1]).intersection(set(df2[seqname2])):
        c1 = df1[cols1][df1[seqname1]==contig].reset_index(drop=True)
        c2 = df2[cols2][df2[seqname2]==contig].sort(start2).drop_duplicates(start2).reset_index(drop=True)
        imerge = pd.cut(c1[end1], c2[start2], labels=False)
        keep = imerge >= 0
        cmerge = pd.concat([c1[keep], c2.iloc[imerge[keep]].reset_index(drop=True)], axis=1).reset_index(drop=True)
        cmerge = cmerge[cmerge[start1] <= cmerge[end2]]
        result = cmerge if result is None else pd.concat([result, cmerge])
    return result


def read_vcf(filename, nrows=None, somatic=False, donor=None, skipcols=None):
    names=['chrom','pos','id','ref','alt','qual','filter','info']
    if somatic:
        names = names + ['format','normal','tumor']
    usecols = None if skipcols==None else [i for i in range(len(names)) if names[i] not in skipcols]
    alldna = pd.read_table(filename, names=names,
                    dtype={'chrom':np.character}, comment='#', index_col=None, nrows=nrows, usecols=usecols)
    if len(alldna.chrom) > 0 and not alldna.iloc[0].chrom.startswith('chr'):
        alldna.chrom = alldna.chrom.apply(lambda c: 'chr'+c)
    if donor != None:
        alldna['donor'] = donor
    return alldna

# find if a position is in a set of regions
# regions must be sorted by start and non-overlapping
def in_regions(pos, regions, start='start', end='end'):
    ix = pd.cut([pos], np.array(zip(regions[start], regions[end] + 1)).reshape(2 * len(regions[start])), labels=False, include_lowest=True, right=False)
    return (ix[0] % 2) == 0
def score_mut(m, mx=None):
    if mx == None:
        mx = jaspar[m['gene']]
    n = len(mx[mx.keys()[0]])
    p = int(m['pos'])
    refseq = hg19[m['chr']][p-n-1:p+n]
    refscore = reduce(max, (score for pos, score in mx.search(refseq, threshold=0.5) if score > 0 and pos >=0 and pos < n), 0)
    mutseq = refseq[:n] + m['var_allele'] + refseq[n+1:]
    mutscore = reduce(max, (score for pos, score in mx.search(mutseq, threshold=0.5) if score > 0 and pos >=0 and pos < n), 0)
    return (mutscore - refscore) / mx.max

def score_mut_all(m):
    max_score = 0
    for g,mx in jaspar.iteritems():
        s = score_mut(m, mx)
        if abs(s) > abs(max_score):
            max_score = s
    return max_score

def maxabs(v):
    return v[np.argmax(np.absolute(v))]

# get enhancer-gene associations

# Ernst et al
def _read_ernst_file(f, fdir):
    df = pd.read_csv(os.path.join(fdir, f), names=['enhchr','enhstart','enhend','gene_name'], sep='\t')
    df['CellLine'] = f.split('_')[1]
    return df
def read_ernst(path='encode/ernst', dir='d:/sequence'):
    return pd.concat((_read_ernst_file(f, os.path.join(dir, path)) for f in os.listdir(os.path.join(dir, path)) if f.startswith('links_') & f.endswith('.txt')))

# Sheffield et al
