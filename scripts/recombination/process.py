#!/bin/python3
#
# Script to run ripples and filtration pipeline on remote GCP machine.

import time
import random
import sys
import subprocess
import datetime
import os
import os.path
import errno

version = sys.argv[1]
mat = sys.argv[2]
start_range = sys.argv[3]
end_range = sys.argv[4]
# Subdir within "results" directory to place results
out = sys.argv[5]
bucket_id = sys.argv[6]
# Remote GCP Storage Bucket location to put end results
results = sys.argv[7]
reference = sys.argv[8]
raw_sequences = sys.argv[9]

pipeline_dir = os.getcwd()

def convert(n):
    return str(datetime.timedelta(seconds = n))

def parse_ripples_command(version, mat, start, end):
    # Expecting ripples output (recombination.txt and descendents.txt)
    # in recombination/filtering to start this pipeline
    command = [version, "-i", mat, "-n", "2", "-S", start, "-E", end, "-d", "filtering/data"]
    return command

# Check starting directory is correct
if (os.path.exists("process.py") == False):
    print()
    print("ERROR: Process.py not found in current directory. Check starting directory in container.") 
    print()
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), "process.py")

# Copy over protobuf, raw sequence file and reference from GCP Storage into irectory inside container
print("Copying input MAT: {} from GCP Storage bucket into local directory on remote machine.".format(mat))
if not os.path.exists(mat):
  subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, mat), pipeline_dir])
if not os.path.exists(reference):
  subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, reference), pipeline_dir])
if not os.path.exists(raw_sequences):
  subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, raw_sequences), pipeline_dir])

# Run ripples on current GCP instance
cmd = parse_ripples_command(version, mat, start_range, end_range)
subprocess.run(cmd)

filtration = ["./run_ripples_filtration.sh", mat, raw_sequences, reference, results, out]
# Run putative recombinants through post-processing filtration pipeline
subprocess.run(filtration)
