#!/usr/bin/python

"""
python script to join inbreeding information and patch age in a file
"""

import csv

def join_coln_age(indfile, patchfile, outfile, header=True):
    """
    """
    csvwriter = csv.writer(outfile)
    ind = csv.reader(indfile)
    patch = csv.reader(patchfile)
    # find the patch columns to keep
    try:
        patch_head = patch.next()
    except StopIteration:
        raise IOError("patch file is empty")
    patch_patch_col = patch_head.index("patch")
    patch_rep_col = patch_head.index("replicate")
    patch_gen_col = patch_head.index("generation")
    patch_age_col = patch_head.index("colnAge")
    # find the individual columns to keep
    try:
        ind_head = ind.next()
    except StopIteration:
        raise IOError("individual file is empty")
    ind_id_col = ind_head.index("id")
    ind_rep_col = ind_head.index("replicate")
    ind_gen_col = ind_head.index("generation")
    ind_patch_col = ind_head.index("patch")
    ind_coa_col = ind_head.index("coancestry")
    ind_fit_col = ind_head.index("fitness")
    # make a patch to age dictionary
    patch2age = {}
    for cols in patch:
        patch_id = cols[patch_patch_col]
        rep = cols[patch_rep_col]
        gen = cols[patch_gen_col]
        age = cols[patch_age_col]
        patch2age[(patch_id, rep, gen)] = age
    # write out the header
    if header:
        row = (
                'replicate',
                'id',
                'patch',
                'generation',
                'coancestry',
                'fitness',
                'colnAge',
                )
        csvwriter.writerow(row)
    # write out the file
    for cols in ind:
        ind_id = cols[ind_id_col]
        rep = cols[ind_rep_col]
        gen = cols[ind_gen_col]
        patch_id = cols[ind_patch_col]
        coa = cols[ind_coa_col]
        fit = cols[ind_fit_col]
        # get the patch age
        age = patch2age[(patch_id, rep, gen)]
        # write to the csv
        row = (
                rep,
                ind_id,
                patch_id,
                gen,
                coa,
                fit,
                age,
                )
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
            outfile = sys.stdout if value == '-' else file(value, 'w')
    if len(arg_ops) != 2:
        sys.stderr.write("usage: python %s [OPTIONS] INDCSV PATCHCSV\n" % \
                (sys.argv[0]))
        sys.exit(1)
    indfile = sys.stdin if arg_ops[0] == '-' else file(arg_ops[0], 'r')
    patchfile = sys.stdin if arg_ops[1] == '-' else file(arg_ops[1], 'r')

    # generate the combined file
    join_coln_age(indfile, patchfile, outfile)
