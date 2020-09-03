#! /usr/bin/env python
import argparse
import os
import json

parser = argparse.ArgumentParser(description="Test JSON parsing")
parser.add_argument('--params')

args = parser.parse_args()

with open(args.params, 'r') as f:
    json_params = json.loads(f.read())

print(f"Contents of {args.params} are: {dir(json_params)}")
