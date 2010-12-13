#!/usr/bin/python

"""
python script to calculate inbreeding depression (of a sort) using a conancestry
and fitness file output by coafitcsv.py
"""

import csv

def inbr_depr_csv(indfile, patchfile, outfile):
    ind_reader = csv.reader(indfile)
    patch_reader = csv.reader(patchfile)
    writer = csv.writer(outfile)
    # get the headers for the individual file
    try:
        headers = ind_reader.next()
    except StopIteration:
        raise IOError("individual CSV file is empty")
    # get the coancestry and fitness column indices
    rep_col = headers.index('replicate')
    gen_col = headers.index('generation')
    patch_col = headers.index('patch')
    coa_col = headers.index('coancestry')
    fit_col = headers.index('fitness')
    # calculate the metapop mean coancestry for each gen and rep
    mean_coa = {}
    count = {}
    for cols in ind_reader:
        rep = cols[rep_col]
        gen = cols[gen_col]
        coa = float(cols[coa_col])
        key = (rep, gen)
        if not key in mean_coa:
            mean_coa[key] = 0
            count[key] = 0
        mean_coa[key] += coa
        count[key] += 1
    for key in mean_coa:
        mean_coa[key] /= count[key]
    # reset the file pointer for rereading
    indfile.seek(0)
    ind_reader.next()
    # get the mean fitness for inbred and outbred classes in each patch
    fitness = {}
    count = {}
    for cols in ind_reader:
        rep = cols[rep_col]
        gen = cols[gen_col]
        patch = cols[patch_col]
        coa = float(cols[coa_col])
        fit = float(cols[fit_col])
        is_inbr = (coa > mean_coa[(rep, gen)])
        key = (rep, gen, patch)
        if not key in fitness:
            fitness[key] = [0, 0]
            count[key] = [0, 0]
        fitness[key][is_inbr] += fit
        count[key][is_inbr] += 1
    # get the header for the patch file
    try:
        headers = patch_reader.next()
    except StopIteration:
        raise IOError("patch CSV file is empty")
    # get the coancestry and fitness column indices
    rep_col = headers.index('replicate')
    gen_col = headers.index('generation')
    patch_col = headers.index('patch')
    # write the headers to the outfile
    row = list(headers) + ['inbrDepr']
    writer.writerow(row)
    # join the difference in mean fitness with each patch and write it out
    for cols in patch_reader:
        rep = cols[rep_col]
        gen = cols[gen_col]
        patch = cols[patch_col]
        key = (rep, gen, patch)
        # get the inbreeding depression
        depr = None
        if key in count:
            if count[key][False] and count[key][True]:
                outbr_fit = fitness[key][False] / count[key][False]
                inbr_fit = fitness[key][True] / count[key][True]
                depr = outbr_fit - inbr_fit
        if depr is None:
            depr = 'NaN'
        # write out the row
        row = list(cols) + [depr]
        writer.writerow(row)

if __name__ == '__main__':
    import sys
    from argproc import ArgProcessor
    # process the arguments
    indfile = sys.stdin
    patchfile = sys.stdin
    outfile = sys.stdout
    flags = {}
    params = {
            'output': ('o', 'output'),
            }
    arg_proc = ArgProcessor(flags, params)
    [flag_ops, param_ops, arg_ops] = arg_proc(sys.argv[1:])
    for param, value in param_ops:
        if param == 'output':
            outfile = sys.stdout if value == '-' else file(value, 'w')
    if len(arg_ops) != 2:
        raise IOError("expected 2 arguments")
    indfile = file(arg_ops[0], 'r')
    patchfile = file(arg_ops[1], 'r')
    # process the file
    inbr_depr_csv(indfile, patchfile, outfile)
