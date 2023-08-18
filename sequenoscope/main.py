#!/usr/bin/env python
import sys


modules = {'analyze': 'map reads to a target and produce a report with sequencing statistics',
            'plot': 'generate plots based on seq manifest files',
            'filter_ONT': 'filter reads from a fastq file based on a sequencing summary file'
            }

module_ordered = ['analyze',
                'plot',
                'filter_ONT'
                ]

def print_usage_and_exit():
    print('Usage: sequenoscope <command> <required arguments>', file=sys.stderr)
    print('\nTo get full help for a command use one of:\nsequenoscope <command> -h\nsequenoscope <command> --help\n', file=sys.stderr)
    print('\nAvailable commands:\n', file=sys.stderr)
    max_module_length = max([len(x) for x in list(modules.keys())]) + 1
    for command in module_ordered:
        print('{{0: <{}}}'.format(max_module_length).format(command), modules[command], sep=' ', file=sys.stderr)
    sys.exit(0)

def main():
    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '-help', '--help']:
        print_usage_and_exit()

    module = sys.argv.pop(1)

    if module not in modules:
        print('Task "' + module + '" not recognised. Cannot continue.\n', file=sys.stderr)
        print_usage_and_exit()

    exec("import sequenoscope.{}.{}".format(module, module))
    exec("sequenoscope.{}.{}".format(module,module) + '.run()')

if __name__ == '__main__':
    main()







