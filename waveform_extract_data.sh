# file name: job.condor
Executable = /home/amedina/CosmicRay.git/trunk/waveform_extract_data.py
output = analysis_data.out
error = analysis_data.err
log = job_data.log
notification = never
request_cpus = 4
request_memory = 3000

#+AccountingGroup="1_week.$ENV(USER)"

#arguments = $(Item)
#Requirements = has_avx =?= true

# use the current metaproject environment
getenv = True

arguments = $(Item)

queue 1 Item from arguments_data.txt

