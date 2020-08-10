#!/bin/bash
source /app/lmod/lmod/init/profile

# get the data
ml Python/3.8.2-GCCcore-9.3.0
ml scanpy

# debugging script.


# call the python script
python convert-to-loom-2.py
