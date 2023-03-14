import numpy as np
import os
import pandas as pd
import config

def generate_alt_codons(codon, loc, alt):
    if loc not in [0, 1, 2]:
        raise ValueError("Location must be 0, 1, or 2")
    prefix = codon[:loc]
    suffix = codon[loc + 1:]

    new_codon = prefix + alt + suffix
    return new_codon


def generate_all_alt_codons(codon, loc):
    basepairs = ['A', 'T', 'C', 'G']
    if loc not in [0, 1, 2]:
        raise ValueError("Location must be 0, 1, or 2")
    prefix = codon[:loc]
    suffix = codon[loc+1:]
    res = []
    for i in range(4):
        new_codon = prefix + basepairs[i] + suffix
        res.append(new_codon)
    return res


def compute_mut_type_dicts():
    basepairs = ['A', 'T', 'C', 'G']
    site_type_dict = {1: '1D', 2:'2D', 3:'3D', 4:'4D'}
    codon_type_dict = {}
    codon_mut_dict = {}
    for codon in codon_dict:
        types = []
        muts = []
        for i in range(3):
            aa = codon_dict[codon]
            alt_codons = generate_all_alt_codons(codon, i)
            alt_type = {}
            for j, alt in enumerate(alt_codons):
                alt_aa = codon_dict[alt]
                if alt_aa==aa:
                    alt_type[basepairs[j]] = 's'
                elif (alt_aa=='*') and (aa!='*'):
                    alt_type[basepairs[j]] = 'n'
                elif (aa=='*') and (alt_aa!='*'):
                    # stop codon got replaced by aa
                    alt_type[basepairs[j]] = 'nn'
                else:
                    alt_type[basepairs[j]] = 'm'
            alt_aas = [aa==codon_dict[alt] for alt in alt_codons]
            num_same = np.sum(alt_aas)
            types.append(site_type_dict[num_same])
            muts.append(alt_type)
        codon_type_dict[codon] = types
        codon_mut_dict[codon] = muts
    return codon_type_dict, codon_mut_dict


def annotate_site_types(seq, strand):
    assert(len(seq)%3 == 0)
    assert(strand in ['+', '-'])
    types = []
    if strand == '-':
        seq = seq.reverse_complement()
    for i in range(len(seq) // 3):
        codon = str(seq[i*3:i*3+3])
        # check if the sequence contain ambiguous nucleotide
        bad_nn = False
        for nn in codon:
            if nn not in allowed_basepairs:
                bad_nn = True
        if bad_nn:
            types.append(['NA', 'NA', 'NA'])
        else:
            types.append(codon_type_dict[codon])
    res = np.hstack(types)
    if strand == '-':
        res = res[::-1]
    return res


allowed_basepairs = {'A', 'T', 'C', 'G'}
codons = pd.read_csv(os.path.join(config.ROOT_DIR, 'dat', 'codons.csv'))
codon_dict = {row['Codon']:row['AA'] for i, row in codons.iterrows()}
codon_type_dict, codon_mut_dict = compute_mut_type_dicts()
