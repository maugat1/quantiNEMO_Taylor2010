"""
Provides the ArgProcessor class.
"""
__author__ = "Malcolm Augat <ma5ke@virginia.edu>"

__all__ = ['ArgProcessor']

class ArgProcessor(object):
    """
    Class for processing file arguments and options.
    """
    def __init__(self, flags={}, params={}):
        """
        Creates an ArgProcessor object for processing the given flag and
        parameter options. flags and params should be dictionaries with option
        keys as ditionary keys, and lists of options corresponding to the key as
        values.
        TODO: explain better and in more depth
        """
        # build a reverse dictionary of the flags
        self.short_flags = {}
        self.long_flags = {}
        for (key, flag_list) in flags.items():
            if key in params:
                raise ValueError(\
                        "option key \"%s\" cannot be both a flag and a param" \
                        % (key))
            for flag in flag_list:
                if len(flag) == 1:
                    if flag in self.short_flags:
                        raise ValueError("repeated option \"%s\"" % (flag))
                    self.short_flags[flag] = key
                else:
                    if flag in self.long_flags:
                        raise ValueError("repeated option \"%s\"" % (flag))
                    self.long_flags[flag] = key
        # build a reverse dictionary of the parameters
        self.short_params = {}
        self.long_params = {}
        for (key, param_list) in params.items():
            for param in param_list:
                if len(param) == 1:
                    if param in self.short_params or param in self.short_flags:
                        raise ValueError("repeated option \"%s\"" % (param))
                    self.short_params[param] = key
                else:
                    if param in self.long_params or param in self.long_flags:
                        raise ValueError("repeated option \"%s\"" % (param))
                    self.long_params[param] = key

    def __call__(self, args):
        """
        Calls the process function.
        """
        return self.process(args)

    def process(self, args):
        """
        For each given argument, attempts to match it against the list of flags
        and parameters used to create the ArgProcessor object and returns
        [flag_ops, param_ops, arguments]. flag_ops is the ordered list of flag
        keys found in args, param_ops is a list of (param key, value) tuples for
        matched param keys and their values, and arguments is a list of
        non-option arguments.

        Example:
        >>> ap = ArgProcessor( flags={'f': ('f', 'flag')}, \
        ... params={'p': ('p', 'param')} )
        >>> [f, p, a] = ap.process( \
        ... ['arg1', '-f', '--param=val1', '-p', 'val2', 'arg2'])
        >>> f
        ['f']
        >>> p
        [('p', 'val1'), ('p', 'val2')]
        >>> a
        ['arg1', 'arg2']
        """
        flags = []
        params = []
        arguments = []
        val_next = None
        for arg in args:
            if not val_next is None:
                # waiting short parameter value
                params.append((self.short_params[val_next], arg))
                val_next = None
            elif not arg:
                # empty argument or option
                raise TypeError("arguments must be non-empty strings")
            elif arg == '-':
                # argument
                arguments.append(arg)
            elif arg[0] == '-':
                # option
                arg = arg[1:]
                if arg[0] == '-':
                    # long option
                    arg = arg[1:]
                    if arg in self.long_flags:
                        # flag
                        flags.append(self.long_flags[arg])
                    else:
                        val = ''
                        if '=' in arg:
                            # option in the form --arg=val
                            [arg, val] = arg.split('=', 1)
                        if arg in self.long_params:
                            # parameter
                            params.append((self.long_params[arg], val))
                        else:
                            # unknown option
                            raise ValueError("unknown option \"--%s\"" % (arg))
                else:
                    # short option
                    for i in range(len(arg)):
                        char = arg[i]
                        if char in self.short_flags:
                            # flag
                            flags.append(self.short_flags[char])
                        elif char in self.short_params:
                            # param
                            if i + 1 == len(arg):
                                val_next = char
                            else:
                                params.append((self.short_params[char], ''))
                        else:
                            # unknown option
                            raise ValueError("unknown option \"-%s\"" % (char))
            else:
                # argument
                arguments.append(arg)
        if val_next:
            params.append((self.short_params[val_next], ''))
        return flags, params, arguments
