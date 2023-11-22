# file name: job.condor
Executable = ./analysis.py
output = analysis.out
error = analysis.err
log = job.log
notification = never
request_cpus = 4
request_memory = 5000

#+AccountingGroup="1_week.$ENV(USER)"

#arguments = $(Item)
#Requirements = has_avx =?= true

# use the current metaproject environment
getenv = True

arguments = $(Item)

queue 1 Item from arguments.txt

