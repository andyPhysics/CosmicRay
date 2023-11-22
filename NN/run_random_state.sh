#!/bin/bash

for i in {0..4..1}
do
    python ./scripts/NN_fitter2.py $i
done
