#!/usr/bin/env python3
"""
Permutation Test for MHA Group Analysis

This script performs permutation tests to calculate empirical p-values for MHA (Microhomology Arms) groups
based on the expected probabilities. The test is conducted for multiple MHA groups (MHA=0, MHA≥1, MHA≥2, MHA≥3),
and p-values are calculated using a binomial distribution with 100,000 simulations.

Usage:
    python permutationTest.py <input_file>

Arguments:
    input_file      : Path to the input file containing TD counts and observed MHA values.

Example:
    python permutationTest.py input_data.txt
"""

import pandas as pd
import numpy as np
from scipy.stats import binom
import sys

def permutation_test(observed_values, n_td_list, expected_p, n_sim=100000, test='greater'):
    """
    Performs a permutation test to calculate the empirical p-value.

    Args:
        observed_values (list): List of observed proportions for each group.
        n_td_list (list): List of TD counts for each species.
        expected_p (float): Expected proportion under the null hypothesis.
        n_sim (int): Number of simulations to perform (default is 100,000).
        test (str): Type of test ('greater' or 'less') to assess deviation from the expected proportion.

    Returns:
        float: The empirical p-value.
    """
    observed_mean = np.mean(observed_values)
    extreme_count = 0
    
    # Perform simulations
    for _ in range(n_sim):
        sim_means = []
        for n_td in n_td_list:
            # Generate binomial random variable for each TD count
            sim_count = binom.rvs(n=n_td, p=expected_p)
            sim_prop = sim_count / n_td
            sim_means.append(sim_prop)
        sim_mean = np.mean(sim_means)
        
        # Count extreme simulations based on test direction
        if test == 'greater' and sim_mean >= observed_mean:
            extreme_count += 1
        elif test == 'less' and sim_mean <= observed_mean:
            extreme_count += 1
    
    # Calculate empirical p-value (with pseudo-counts)
    p_value = (extreme_count + 1) / (n_sim + 1)
    return p_value

def main(input_file):
    """
    Main function to perform permutation tests for MHA groups and print results.
    
    Args:
        input_file (str): Path to the input data file.
    """
    # Load the input data
    df = pd.read_csv(input_file, sep='\t', header=None)
    n_td_list = df[0].tolist()  # First column: TD counts
    
    # Define MHA groups parameters (column index, expected probability, test direction)
    mha_groups = [
        (1, 0.75, 'less'),     # MHA=0bp vs 75% (expect lower)
        (2, 0.25, 'greater'),  # MHA≥1bp vs 25% (expect higher)
        (3, 0.0625, 'greater'),# MHA≥2bp vs 6.25% 
        (4, 0.015625, 'greater') # MHA≥3bp vs 1.56%
    ]
    
    results = []
    
    # Perform permutation tests for each MHA group
    for col_idx, expected_p, test_dir in mha_groups:
        # Convert observed percentages to decimal
        obs_values = df[col_idx].str.rstrip('%').astype(float).div(100).tolist()
        
        # Perform permutation test
        p_val = permutation_test(
            observed_values=obs_values,
            n_td_list=n_td_list,
            expected_p=expected_p,
            test=test_dir
        )
        
        # Store the results
        results.append({
            'MHA Group': ['0bp', '≥1bp', '≥2bp', '≥3bp'][col_idx-1],
            'Expected Value': f"{expected_p*100:.2f}%",
            'Test Direction': test_dir,
            'Empirical p-value': f"{p_val:.4e}"
        })
    
    # Print the results in a tabular format
    print("\n{:<8} {:<12} {:<14} {:<15}".format(
        "Group", "Expected", "Test Direction", "p-value"))
    print("-" * 45)
    for res in results:
        print("{:<8} {:<12} {:<14} {:<15}".format(
            res['MHA Group'],
            res['Expected Value'],
            res['Test Direction'],
            res['Empirical p-value']
        ))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python permutationTest.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    main(input_file)

