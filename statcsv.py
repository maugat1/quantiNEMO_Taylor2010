#!/usr/bin/python

"""
python script to cleanup the CSV stat files created by quantiNEMO
"""

import csv
import re

header_p = re.compile(r'^(.*?)(_p(\d+))?(_t\d+)?$', re.I)

def stat_csv(statfile, outfile):
    """
    """
    stat = csv.reader(statfile)
    out = csv.writer(outfile)
    # read the header
    try:
        headers = stat.next()
    except StopIteration:
        raise IOError("stat CSV file is empty")
    # build maps of columns to save
    i2header = {}
    header2col = {}
    col2header = {}
    next_col = 0
    i2patch = {}
    patch2i = {}
    for i in range(len(headers)):
        m = header_p.match(headers[i])
        header = m.groups()[0]
        trait_num = m.groups()[3]
        if not trait_num is None:
            header += trait_num
        i2header[i] = header
        patch_num = m.groups()[2]
        if not patch_num is None:
            patch_num = int(patch_num)
            i2patch[i] = patch_num
            if not patch_num in patch2i:
                patch2i[patch_num] = []
            patch2i[patch_num].append(i)
        if not header in header2col:
            header2col[header] = next_col
            col2header[next_col] = header
            next_col += 1
    # write out the header
    if len(i2patch):
        header2col['patch'] = next_col
        col2header[next_col] = 'patch'
        next_col += 1
    headers = []
    for i in range(len(col2header)):
        headers.append(col2header.get(i, ''))
    out.writerow(headers)
    # process the file
    for cols in stat:
        # get the shared (non-patch) values
        base_row = [''] * len(headers)
        for i in range(len(cols)):
            if not i in i2patch:
                base_row[header2col[i2header[i]]] = cols[i]
        # check if there are no patch values
        if not len(i2patch):
            out.writerow(base_row)
            continue
        # get the patch specific stuff
        rows = {}
        for i in range(len(cols)):
            if i in i2patch:
                patch_num = i2patch[i]
                if not patch_num in rows:
                    rows[patch_num] = [x for x in base_row]
                    rows[patch_num][header2col['patch']] = patch_num
                rows[patch_num][header2col[i2header[i]]] = cols[i]
        for row in rows.values():
            out.writerow(row)

if __name__ == '__main__':
    import sys
    from argproc import ArgProcessor
    # process the arguments
    statfile = sys.stdin
    outfile = sys.stdout
    flags = {}
    params = {}
    arg_proc = ArgProcessor(flags, params)
    [flag_ops, param_ops, arg_ops] = arg_proc(sys.argv[1:])
    if len(arg_ops) != 2:
        raise IOError("expected 2 arguments")
    statfile = sys.stdin if arg_ops[0] == '-' else file(arg_ops[0], 'r')
    outfile = sys.stdout if arg_ops[1] == '-' else file(arg_ops[1], 'w')
    # process the files
    stat_csv(statfile, outfile)
