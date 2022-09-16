import sys
import os
import subprocess
import json
from typing import *


def slice_project_args(args: List[str] = sys.argv) -> List[str]:
    if '--run' in args:
        args = args[args.index('--run') + 2:]
    else:
        args = args[1:]
    return list(filter(lambda x: x != "-v" and x != '--verbose', args))


def run_one(run: str, verbose: bool = False, **kwargs):
    dir_base = os.path.abspath(os.getcwd())
    dir_run = os.path.abspath(run)
    os.chdir(dir_run)
    os.environ['ARGS'] = json.dumps(sys.argv)
    if verbose:
        print(f"os.environ['ARGS'] = {os.environ['ARGS']}")
    subprocess.call(['jupyter-nbconvert', '--to',
                    'html', '--execute', f"{run}.ipynb"])
    os.chdir(dir_base)


def done(**kwargs):
    print("=" * 20, "DONE", "=" * 20)
    if isinstance(kwargs, dict) and len(kwargs.keys()) > 0:
        print("RESULT:", kwargs)
