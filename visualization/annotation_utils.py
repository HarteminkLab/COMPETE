
'''
utils related to visualize annotations
'''

import roman
import pandas as pd
import sys

def find_name(x):
    info = dict([_.split('=') for _ in x.split(';')])
    return info['gene'] if 'gene' in info else info['Name']

def get_orf_data(chromo, start, end, orf_annotation = None):
    
    if orf_annotation is None:
        sys.exit("For ORF annotation, a GFF file needs to be provided. Such as this one: https://github.com/jianlingzhong/COMPETE_examples/blob/master/saccharomyces_cerevisiae.20080621.gff")

    genes = pd.read_csv(orf_annotation, sep = '\t', comment='#', header = None)

    # only need the gene features in the gff
    genes = genes.loc[genes.iloc[:, 2] == 'gene', ]

    # don't need Mito and 2-micron chromosome features
    genes = genes.loc[(genes.iloc[:, 0] != 'chrMito') & (genes.iloc[:, 0] != '2-micron')]

    # just need the following columns
    genes = genes.iloc[:, [0, 3,4,6,8]]
    genes.rename(columns={0:'chromo', 3:'start', 4:'end', 6:'strand', 8:'name'}, inplace = True)

    # convert 'chrI' to 1
    genes.chromo = genes.chromo.apply(lambda x: roman.fromRoman(x[3:]))

    genes = genes.loc[(genes.chromo == chromo)]
    selection = ~((genes.end < start) | (genes.start > end))
    genes = genes.loc[selection]
    genes.name = genes.name.apply(find_name)

    return genes

def get_macisaac_data(chromo, start, end, macisaac_file = None):
    '''read the macisaac binding sites data. currently only data with p-value 0.05 and conservation level 2 is included'''

    if macisaac_file is None:
        sys.exit("For Macisaac binding sites annotation, a supplemental file from Macisaac et al. needs to be provided. Such as this one: https://github.com/jianlingzhong/COMPETE_examples/blob/master/macisaac_p005_c2.csv")

    binding = pd.read_csv(macisaac_file, sep = '\t')
    binding = binding.loc[binding.chromo == chromo]

    selection = ~((binding.end < start) | (binding.start > end))
    binding = binding.loc[selection]

    return binding
