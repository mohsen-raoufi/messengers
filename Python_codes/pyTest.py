"""
Title: Runs a Single Simulation to Test

Description: This script runs a test simulation to be used for testing the implementation.

Author: Mohsen Raoufi

Contact: mohsenraoufi@icloud.com

Affiliation: Research Cluster of Excellence, "Science of Intelligence"

Date Created: August, 2024

Version: 1.0

Usage: Run this script with Python 3.8 or later. Ensure you create a "data" directory in the working directory before execution.

License: Distributed under the MIT License. See LICENSE.txt for more information.

Citation: Please cite our work if you use this code:
"Messengers: Breaking Echo Chambers in Collective Opinion Dynamics with Homophily"
"""

import run_one_config as roc

print("Running one experiment from run_one_config.py", flush=True)
result, params = roc.__main__()
print("finished running one experiment")