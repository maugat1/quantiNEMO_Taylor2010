#!/usr/bin/python

"""
python script to calculate parental coancestry of individuals
"""

import csv
import re

def coancestry(loci1, loci2):
    """
    Arguments:
        ind1: list of allele index tuples
        ind2: list of allele index tuples
    Returns:
        coancestry: a float representing the degree of similarity between the
            two individuals
    """
    # get and check the number of loci
    nb_loci = len(loci1);
    if nb_loci != len(loci2):
        raise IndexError("loci sets are not the same length: %d and %d" % \
                (nb_loci, len(loci2)))
    # calculate the number of matching alleles
    coa = 0
    for i in range(nb_loci):
        loc1 = loci1[i]
        loc2 = loci2[i]
        coa += (loc1[0] == loc2[0]) + \
                (loc1[0] == loc2[1]) + \
                (loc1[1] == loc2[0]) + \
                (loc1[1] == loc2[1])
    # divide by the maximum possible number of matches
    coa /= (4.0 * nb_loci)
    return coa

locus_p = re.compile(r'(trait-\d\d*)?_?loc(us)?-\d\d*', re.I)

def coa_fit_csv(phenfile, markfile, outfile, header=True):
    """
    """
    csvwriter = csv.writer(outfile)
    phen = csv.reader(phenfile)
    mark = csv.reader(markfile)
    # find the phenotype columns to keep
    try:
        phen_head = phen.next()
    except StopIteration:
        raise IOError("phenotype file is empty")
    phen_id_col = phen_head.index("id")
    phen_patch_col = phen_head.index("patch")
    phen_mom_col = phen_head.index("mother_id")
    phen_dad_col = phen_head.index("father_id")
    phen_fit_col = phen_head.index('fitness')
    try:
        phen_rep_col = phen_head.index('replicate')
    except IndexError:
        phen_rep_col = None
    # find the marker columns to keep
    try:
        mark_head = mark.next()
    except StopIteration:
        raise IOError("marker file is empty")
    mark_id_col = mark_head.index("id")
    mark_loc_cols = []
    for i in range(len(mark_head)):
        if locus_p.match(mark_head[i]):
            mark_loc_cols.append(i)
    try:
        mark_rep_col = mark_head.index('replicate')
    except IndexError:
        mark_rep_col = None
    # make a parent genotype dictionary from marker file
    parent_loci = {}
    for cols in mark:
        ind = cols[mark_id_col]
        loci = []
        for i in mark_loc_cols:
            a1, a2 = cols[i].split(':')
            loci.append((int(a1), int(a2)))
        rep = 1 if mark_rep_col is None else cols[mark_rep_col]
        parent_loci[(ind, rep)] = loci
    # write out the header
    if header:
        row = ('replicate', 'id', 'patch', 'coancestry', 'fitness')
        csvwriter.writerow(row)
    # write out coancestry and fitness
    for cols in phen:
        ind = cols[phen_id_col]
        patch = cols[phen_patch_col]
        rep = 1 if phen_rep_col is None else cols[phen_rep_col]
        mom = (cols[phen_mom_col], rep)
        dad = (cols[phen_dad_col], rep)
        fit = cols[phen_fit_col]
        # get the parents' coancestry
        coa = coancestry(parent_loci[mom], parent_loci[dad])
        # write to the csv
        row = (rep, ind, patch, coa, fit)
        csvwriter.writerow(row)

if __name__ == '__main__':
    import sys
    from argproc import ArgProcessor

    # process the arguments
    outfile = sys.stdout

    flags = {}
    params = {
            'output': ('o', 'output'),
            }
    arg_proc = ArgProcessor(flags, params)
    [flag_ops, param_ops, arg_ops] = arg_proc(sys.argv[1:])
    for param, value in param_ops:
        if param == 'output':
            outfile = file(value, 'w')
    if len(arg_ops) < 2:
        sys.stderr.write("usage: python %s [OPTIONS] PHENCSV MARKCSV\n")
        sys.exit(1)
    phenfile = file(arg_ops[0], 'r')
    markfile = file(arg_ops[1], 'r')

    # generate the combined file
    coa_fit_csv(phenfile, markfile, outfile)
