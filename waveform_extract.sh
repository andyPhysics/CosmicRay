# file name: job.condor
Executable = /home/amedina/CosmicRay.git/trunk/waveform_extract_2.py
output = analysis_new.out
error = analysis_new.err
log = job_new.log
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

