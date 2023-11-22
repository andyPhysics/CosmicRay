# file name: job.condor
Executable = ./analysis_data.py
output = analysis_data2.out
error = analysis_data2.err
log = job_data2.log
notification = never
request_cpus = 4
request_memory = 6000

#+AccountingGroup="1_week.$ENV(USER)"

#arguments = $(Item)
#Requirements = has_avx =?= true

# use the current metaproject environment
getenv = True

arguments = $(Item)

queue 1 Item from files.txt

