#!/usr/bin/python

"""
python script to extract information from the FSTAT files output by quantiNEMO
and to rewrite them as csv files.
"""

import csv

def fstat2csv(fstatfile, csvfile, types=None, patch=True, age=True, sex=True, \
        ind_id=True, mom_id=True, dad_id=True, fit=True, default="", \
        header=True, extra=[]):
    """
    Reads an extended FSTAT file from quantiNEMO and uses it to write a CSV file
    with similar information.

    Arguments:
        fstatfile: file reading an extended FSTAT file from quantiNEMO
        csvfile: file writing a CSV file
        types: list of numbers for which types to keep in the CSV (default all)
        patch: whether to keep the individual's patch in the CSV (default yes)
        age: whether to keep the individual's age in the CSV (default yes)
        sex: whether to keep the individual's sex in the CSV (default yes)
        ind_id: whether to keep the individual's id in the CSV (default yes)
        mom_id: whether to keep the individual's mom_id in the CSV (default yes)
        dad_id: whether to keep the individual's dad_id in the CSV (default yes)
        fit: whether to keep the individual's fitness in the CSV (default yes)
        default: default value for when a column is missing
        header: whether to print a header row in the CSV (default yes)
        extra: (name_str, value) tuples of extra columns to add to the CSV file
    """
    # read the first line to get the size of things
    try:
        line = fstatfile.next()
    except StopIteration:
        raise IOError("FSTAT file is empty")
    info = [int(s.strip()) for s in line.split()]
    if len(info) == 2:
        nb_patches, nb_types = info
    elif len(info) == 4:
        nb_patches, nb_types, max_subtype, len_subtype = info
    else:
        raise ValueError("expected 2 or 4 info columns, got %d" % (len(info)))
    # create the csv writer
    csvwriter = csv.writer(csvfile)
    # read the patch names
    if types is None:
        types = range(1, nb_types+1)
    if len(types):
        type_names = dict(zip(types, ("type%d" % (i) for i in types)))
    for i in range(1, nb_types+1):
        try:
            line = fstatfile.next()
        except StopIteration:
            raise IOError("unexpected FSTAT file termination in " + \
                    "type names")
        if len(types):
            type_names[i] = line.strip()
    # write out the header row and get the column indices to keep
    row = []
    cols = []
    if patch:
        row.append("patch")
        cols.append(0)
    if len(types):
        for i in types:
            row.append(type_names[i])
            cols.append(i)
    if age:
        row.append("age")
        cols.append(nb_types + 1)
    if sex:
        row.append("sex")
        cols.append(nb_types + 2)
    if ind_id:
        row.append("id")
        cols.append(nb_types + 3)
    if mom_id:
        row.append("mother_id")
        cols.append(nb_types + 4)
    if dad_id:
        row.append("father_id")
        cols.append(nb_types + 5)
    if fit:
        row.append("fitness")
        cols.append(nb_types + 6)
    for name, value in extra:
        row.append(name)
    if header:
        csvwriter.writerow(row)
    # read the values for each individual
    for line in fstatfile:
        values = line.split()
        row = []
        # append the values we are keeping
        for i in cols:
            if len(info) == 4 and i >= 1 and i <= nb_types:
                values[i] = values[i][:len_subtype] + ':' + \
                        values[i][len_subtype:]
            try:
                row.append(values[i])
            except IndexError:
                row.append(default)
        # append the extra values
        for name, value in extra:
            row.append(value)
        # write out the row to the file
        csvwriter.writerow(row)

if __name__ == '__main__':
    import sys
    from argproc import ArgProcessor
    # process the arguments
    types = None
    patch = age = sex = ind_id = mom_id = dad_id = fit = header = True
    default = ""
    extra = []
    fstatfiles = [sys.stdin]
    csvfile = sys.stdout
    flags = {
            'patch': ('p', 'patch'),
            'nopatch': ('P', 'no-patch'),
            'age': ('a', 'age'),
            'noage': ('A', 'no-age'),
            'sex': ('s', 'sex'),
            'nosex': ('S', 'no-sex'),
            'individual': ('i', 'individual'),
            'noindividual': ('I', 'no-individual'),
            'mother': ('m', 'mother'),
            'nomother': ('M', 'no-mother'),
            'father': ('f', 'father'),
            'nofather': ('F', 'no-father'),
            'fitness': ('w', 'fitness'),
            'nofitness': ('W', 'no-fitness'),
            'header': ('h', 'header'),
            'noheader': ('H', 'no-header'),
            }
    params = {
            'output': ('o', 'output'),
            'types': ('z', 'types'),
            'default': ('d', 'default'),
            'extra': ('x', 'extra'),
            }
    arg_proc = ArgProcessor(flags, params)
    [flag_ops, param_ops, arg_ops] = arg_proc(sys.argv[1:])
    for flag in flag_ops:
        if flag == 'patch':
            patch = True
        elif flag == 'nopatch':
            patch = False
        elif flag == 'age':
            age = True
        elif flag == 'noage':
            age = False
        elif flag == 'sex':
            sex = True
        elif flag == 'nosex':
            sex = False
        elif flag == 'individual':
            ind_id = True
        elif flag == 'noindividual':
            ind_id = False
        elif flag == 'mother':
            mom_id = True
        elif flag == 'nomother':
            mom_id = False
        elif flag == 'father':
            dad_id = True
        elif flag == 'nofather':
            dad_id = False
        elif flag == 'fitness':
            fit = True
        elif flag == 'nofitness':
            fit = False
        elif flag == 'header':
            header = True
        elif flag == 'noheader':
            header = False
    for (param, value) in param_ops:
        if param == 'output':
            csvfile = sys.stdout if value == '-' else file(value, 'w')
        elif param == 'types':
            types = []
            if value.strip():
                types = [int(s.strip()) for s in value.split(',')]
        elif param == 'default':
            default = value
        elif param == 'extra':
            extra = []
            if value.strip():
                for pair in value.split(','):
                    name, val = pair.split(':')
                    extra.append((name.strip(), val.strip()))
    if len(arg_ops):
        fstatfiles = []
    for arg in arg_ops:
        fstatfiles.append(sys.stdin if arg == '-' else file(arg, 'r'))
    # process the files
    for fstatfile in fstatfiles:
        fstat2csv(fstatfile=fstatfile, csvfile=csvfile, types=types, \
                patch=patch, age=age, sex=sex, ind_id=ind_id, mom_id=mom_id, \
                dad_id=dad_id, fit=fit, default=default, header=header, \
                extra=extra)
        header = False
