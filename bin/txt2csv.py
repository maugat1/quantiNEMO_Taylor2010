#!/usr/bin/python

"""
python script to write a CSV file from a tab separated text file created by
quantiNEMO
"""

import csv

def txt2csv(txtfile, csvfile):
    writer = csv.writer(csvfile)
    for line in txtfile:
        writer.writerow(line.strip().split())

if __name__ == '__main__':
    import sys
    from argproc import ArgProcessor
    # process the arguments
    txtfile = sys.stdin
    outfile = sys.stdout
    flags = {}
    params = {}
    arg_proc = ArgProcessor(flags, params)
    [flag_ops, param_ops, arg_ops] = arg_proc(sys.argv[1:])
    if len(arg_ops) != 2:
        raise IOError("expected 2 arguments")
    txtfile = file(arg_ops[0], 'r')
    outfile = file(arg_ops[1], 'w')
    # process the files
    txt2csv(txtfile, outfile)
