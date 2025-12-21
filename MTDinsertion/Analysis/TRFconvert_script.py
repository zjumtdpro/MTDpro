#!/usr/bin/env python3
"""
TRF Output Processor
Version: 1.2
Function: Processes TRF output files to standardized format
"""

import argparse
import glob
import os
from pathlib import Path

def process_trf(input_pattern: str, output_suffix: str = "-tr.dat") -> None:
    """Main processing function"""
    for input_file in glob.glob(input_pattern):
        output_path = Path(input_file).with_suffix(output_suffix)
        
        with open(input_file, 'r') as fin, open(output_path, 'w') as fout:
            current_seq = ""
            for line in fin:
                if line.startswith("Sequence:"):
                    current_seq = line.split()[1].strip()
                elif line[0].isdigit():
                    fout.write(f"{current_seq}\t{line}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process TRF output files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", default="*.dat",
                        help="Input file pattern")
    parser.add_argument("-o", "--output-suffix", default="-tr.dat",
                        help="Output file suffix")
    
    args = parser.parse_args()
    
    process_trf(
        input_pattern=args.input,
        output_suffix=args.output_suffix
    )
