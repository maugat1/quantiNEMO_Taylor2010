#!/usr/bin/python

import sys

def logtime(start, stop, step, ostream=sys.stdout):
    never = stop + 1
    ostream.write("(1 %d" % (never))
    for gen in range(start, never, step):
        if gen != 1:
            ostream.write(", %d 1, %d %d" % (gen - 1, gen, never))
    ostream.write(")\n")

if __name__ == '__main__':
    from argproc import ArgProcessor

    # process the arguments
    ostream = sys.stdout

    flags = {}
    params = {
            'output': ('o', 'output'),
            }
    arg_proc = ArgProcessor(flags, params)
    [flag_ops, param_ops, arg_ops] = arg_proc(sys.argv[1:])
    for param, value in param_ops:
        if param == 'output':
            ostream = sys.stdout if value == '-' else file(value, 'w')
    if len(arg_ops) != 3:
        sys.stderr.write("usage: python %s [OPTIONS] START STOP STEP\n" % \
                (sys.argv[0]))
        sys.exit(1)
    start = int(arg_ops[0])
    stop = int(arg_ops[1])
    step = int(arg_ops[2])

    # write the log times
    logtime(start, stop, step, ostream)

