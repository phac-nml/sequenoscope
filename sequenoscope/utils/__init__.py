#!/usr/bin/env python
from subprocess import Popen, PIPE
import os
import hashlib

def run_command(command):
    p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')
    return stdout, stderr

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def compute_sha256(file_name):
    hash_sha256 = hashlib.sha256()
    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(5000), b""):
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()  

def format_time(seconds):
    """Returns time in the most appropriate unit."""
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    elif seconds < 3600:
        return f"{seconds / 60:.2f} minutes"
    else:
        return f"{seconds / 3600:.2f} hours"