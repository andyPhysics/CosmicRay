# file name: job.condor
Executable = ./scripts/extract_data.py
output = extract_data.out
error = extract_data.err
log = job.log
notification = never
request_cpus = 4
request_memory = 10000

#arguments = $(Item)
#Requirements = has_avx =?= true

# use the current metaproject environment
getenv = True

arguments = $(Item)

queue 1 Item from arguments.txt

