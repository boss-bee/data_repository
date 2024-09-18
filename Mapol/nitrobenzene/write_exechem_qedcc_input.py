import json
import sys

input_options_fn = sys.argv[1]

# read input options
with open("input_options.json", 'r') as f:
    input_options = json.load(f)

print(input_options)


