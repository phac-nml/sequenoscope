#!/usr/bin/env python
import sys
import subprocess
import os
import argparse
from sequenoscope.version import __version__
from sequenoscope.utils.__init__ import format_time

modules = {
    'analyze': 'map reads to a target and produce a report with sequencing statistics',
    'plot': 'generate plots based on seq manifest files',
    'filter_ONT': 'filter reads from a fastq file based on a sequencing summary file'
}

module_ordered = ['analyze', 'plot', 'filter_ONT']

def print_usage_and_exit():
    print('Usage: sequenoscope <command> <required arguments>', file=sys.stderr)
    print('\nTo get full help for a command use one of:\nsequenoscope <command> -h\nsequenoscope <command> --help\n', file=sys.stderr)
    print('\nAvailable commands:\n', file=sys.stderr)
    max_module_length = max([len(x) for x in list(modules.keys())]) + 1
    for command in module_ordered:
        print('{{0: <{}}}'.format(max_module_length).format(command), modules[command], sep=' ', file=sys.stderr)
    print('\nOther options:\n', file=sys.stderr)
    print('--check_dependencies  Check if external dependencies (fastp, minimap2, samtools, mash, seqtk) and required Python packages (pysam, plotly) are available', file=sys.stderr)
    print('-v, --version         Show the version and exit', file=sys.stderr)
    print('-h, --help            Show this help message and exit', file=sys.stderr)
    sys.exit(0)

def check_dependencies():
    """
    Check if external dependencies and required Python packages are available.
    """
    # External command-line tools
    required_tools = ["fastp", "minimap2", "samtools", "mash", "seqtk"]

    # Python packages
    required_packages = ["pysam", "plotly"]

    missing_tools = []
    for tool in required_tools:
        if not is_tool_available(tool):
            missing_tools.append(tool)

    missing_packages = []
    for pkg in required_packages:
        if not is_package_available(pkg):
            missing_packages.append(pkg)

    # Report missing tools/packages if any
    if missing_tools or missing_packages:
        if missing_tools:
            print("ERROR: The following required external tools are not found in your PATH:", file=sys.stderr)
            for t in missing_tools:
                print(f" - {t}", file=sys.stderr)
            print("Please install these tools before running Sequenoscope.\n", file=sys.stderr)

        if missing_packages:
            print("ERROR: The following required Python packages are not installed:", file=sys.stderr)
            for p in missing_packages:
                print(f" - {p}", file=sys.stderr)
            print("Please install these Python packages (e.g., via pip) before running Sequenoscope.", file=sys.stderr)

        sys.exit(1)
    else:
        print("All required external tools and Python packages are available.", file=sys.stderr)
        sys.exit(0)

def is_tool_available(tool_name):
    """
    Check if a tool is available in the system PATH by attempting to run 'which <tool>'.
    Returns True if found, False otherwise.
    """
    try:
        subprocess.run(["which", tool_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except subprocess.CalledProcessError:
        return False

def is_package_available(package_name):
    """
    Check if a Python package is available by trying to import it.
    Returns True if import succeeds, False otherwise.
    """
    try:
        __import__(package_name)
        return True
    except ImportError:
        return False

def main():
    # If no arguments or asking for help
    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
        print_usage_and_exit()

    # Version check
    if len(sys.argv) == 2 and sys.argv[1] in ['-v', '--version']:
        print(__version__)
        sys.exit(0)

    # Check dependencies only if requested
    if len(sys.argv) == 2 and sys.argv[1] == '--check_dependencies':
        check_dependencies()

    # Otherwise assume the first argument is a module command
    module = sys.argv.pop(1)

    if module not in modules:
        print(f'Task "{module}" not recognised. Cannot continue.\n', file=sys.stderr)
        print_usage_and_exit()

    exec("import sequenoscope.{}.{}".format(module, module))
    exec("sequenoscope.{}.{}".format(module, module) + '.run()')

if __name__ == '__main__':
    main()
