#!/usr/bin/python

"""
python script to calculate parental coancestry of individuals
"""

import tempfile
import shutil
import re
import csv
import buzhug

def coancestry(loci1, loci2):
    """
    Arguments:
        loci1: tuple of two allele strings
        loci2: tuple of two allele strings
    Returns:
        coancestry: a float representing the degree of similarity between the
            two individuals
    """
    # get and check the number of loci
    if len(loci1) != 2 or len(loci2) != 2:
        raise IndexError("expected two sets of alleles for each parent")
    nb_loci = len(loci1[0]);
    if nb_loci != len(loci1[1]):
        raise IndexError("loci sets are not the same length: %d and %d" % \
                (nb_loci, len(loci1[1])))
    if nb_loci != len(loci1[0]):
        raise IndexError("loci sets are not the same length: %d and %d" % \
                (nb_loci, len(loci2[0])))
    if nb_loci != len(loci2[1]):
        raise IndexError("loci sets are not the same length: %d and %d" % \
                (nb_loci, len(loci2[1])))
    # calculate the number of matching alleles
    coa = 0
    for i in range(nb_loci):
        coa += (loci1[0][i] == loci2[0][i]) + \
                (loci1[0][i] == loci2[1][i]) + \
                (loci1[1][i] == loci2[0][i]) + \
                (loci1[1][i] == loci2[1][i])
    # divide by the maximum possible number of matches
    coa /= (4.0 * nb_loci)
    return coa

locus_p = re.compile(r'(trait-\d\d*)?_?loc(us)?-\d\d*', re.I)

def coa_fit_csv(phenfile, markfile, outfile, dbpath=None, header=True):
    """
    Creates a csv with coancestry and fitness information from the given
    phenotype and marker CSV files. In order to link the two, a database is
    created in the given dbpath (by default a temporary directory is created).
    """
    csvwriter = csv.writer(outfile)
    phen = csv.reader(phenfile)
    mark = csv.reader(markfile)
    # set up the marker database
    if dbpath is None:
        tmppath = tempfile.mkdtemp()
        shutil.rmtree(tmppath, ignore_errors=True)
        parent_db = buzhug.Base(tmppath)
    else:
        parent_db = buzhug.Base(dbpath)
    # find the phenotype columns to keep
    try:
        phen_head = phen.next()
    except StopIteration:
        raise IOError("phenotype file is empty")
    phen_id_col = phen_head.index("id")
    phen_patch_col = phen_head.index("patch")
    phen_mom_col = phen_head.index("mother_id")
    phen_dad_col = phen_head.index("father_id")
    phen_fit_col = phen_head.index("fitness")
    try:
        phen_rep_col = phen_head.index("replicate")
    except IndexError:
        phen_rep_col = None
    try:
        phen_gen_col = phen_head.index("generation")
    except IndexError:
        phen_gen_col = None
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
        mark_rep_col = mark_head.index("replicate")
    except IndexError:
        mark_rep_col = None
    # make a parent genotype database from the marker file
    parent_db.create(("id", str), ("replicate", str, '1'), ("alleles1", str),
            ("alleles2", str), mode="override")
    for cols in mark:
        ind = cols[mark_id_col]
        loci1 = []
        loci2 = []
        for i in mark_loc_cols:
            a1, a2 = cols[i].split(':')
            loci1.append(int(a1) - 1)
            loci2.append(int(a2) - 1)
        alleles1 = ''.join((chr(l1) for l1 in loci1))
        alleles2 = ''.join((chr(l2) for l2 in loci2))
        if mark_rep_col is None:
            parent_db.insert(id=ind, alleles1=alleles1, alleles2=alleles2)
        else:
            parent_db.insert(id=ind, alleles1=alleles1, alleles2=alleles2,
                    replicate=cols[mark_rep_col])
    # write out the header
    if header:
        row = (
                'replicate',
                'generation',
                'id',
                'patch',
                'coancestry',
                'fitness',
                )
        csvwriter.writerow(row)
    # write out coancestry and fitness
    for cols in phen:
        ind = cols[phen_id_col]
        patch = cols[phen_patch_col]
        rep = 1 if phen_rep_col is None else cols[phen_rep_col]
        gen = 0 if phen_gen_col is None else cols[phen_gen_col]
        mom_id = cols[phen_mom_col]
        dad_id = cols[phen_dad_col]
        fit = cols[phen_fit_col]
        # get the parents' coancestry
        mom_record = parent_db.select(["alleles1", "alleles2"], id=mom_id)[0]
        dad_record = parent_db.select(["alleles1", "alleles2"], id=dad_id)[0]
        coa = coancestry((mom_record.alleles1, mom_record.alleles2),
                (dad_record.alleles1, dad_record.alleles2))
        # write to the csv
        row = (
                rep,
                gen,
                ind,
                patch,
                coa,
                fit
                )
        csvwriter.writerow(row)
    # delete the temporary database if necessary
    if dbpath is None:
        shutil.rmtree(tmppath, ignore_errors=True)

if __name__ == '__main__':
    import sys
    from argproc import ArgProcessor

    # process the arguments
    outfile = sys.stdout
    dbpath = None

    flags = {}
    params = {
            'output': ('o', 'output'),
            'dbpath': ('d', 'dbpath'),
            }
    arg_proc = ArgProcessor(flags, params)
    [flag_ops, param_ops, arg_ops] = arg_proc(sys.argv[1:])
    for param, value in param_ops:
        if param == 'output':
            outfile = sys.stdout if value == '-' else file(value, 'w')
        elif param == 'dbpath':
            dbpath = value
    if len(arg_ops) != 2:
        sys.stderr.write("usage: python %s [OPTIONS] PHENCSV MARKCSV\n" % \
                (sys.argv[0]))
        sys.exit(1)
    phenfile = sys.stdin if arg_ops[0] == '-' else file(arg_ops[0], 'r')
    markfile = sys.stdin if arg_ops[1] == '-' else file(arg_ops[1], 'r')

    # generate the combined file
    coa_fit_csv(phenfile, markfile, outfile, dbpath=dbpath)
