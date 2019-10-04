# file name: job.condor
Executable = /home/amedina/CosmicRay/NN/NN.py
output = output.out
error = output.err
log = job.log
notification = never

# use the current metaproject environment
getenv = True

Requirements = has_avx =?= true

queue
