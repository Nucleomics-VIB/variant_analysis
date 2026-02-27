#!/usr/bin/env python3

# script: parse_err.py
# parse clara parabricks .err file and extract the description and times of each subtask
# report as csv for plotting
# SP@NC; 2025-01-21; v1.0

import sys
import re

def extract_info(block):
    program = re.search(r'Program:\s+(.*?)\s+\|\|', block)
    total_time = re.search(r'Total Time:\s+(.*?)\s+\|\|', block)

    if program and total_time:
        program_name = program.group(1).strip()
        time_str = total_time.group(1).strip()
        total_seconds = convert_to_seconds(time_str)
        return program_name, time_str, total_seconds
    return None, None, None

def convert_to_seconds(time_str):
    parts = time_str.split()
    total_seconds = 0
    for i in range(0, len(parts), 2):
        value = int(parts[i])
        unit = parts[i+1]
        if unit.startswith('minute'):
            total_seconds += value * 60
        elif unit.startswith('second'):
            total_seconds += value
    return total_seconds

def process_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
        blocks = content.split('------------------------------------------------------------------------------')

        for block in blocks:
            program, time, seconds = extract_info(block)
            if program and time:
                print(f'"{program}", "{time}", {seconds}')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    process_file(input_file)

