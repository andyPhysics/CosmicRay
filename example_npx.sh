# file name: job.condor
Executable = /home/amedina/CosmicRay.git/trunk/curvature.py
output = curvature.out
error = curvature.err
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

