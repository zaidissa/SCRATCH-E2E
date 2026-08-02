#!/usr/bin/env python3

####!/opt/conda/bin/python

import os
import sys
import re

def is_already_10x_format(sample_name, filename):
    """
    Determine if the filename matches the expected 10x Genomics format.
    {sample_name}_S\d+_L\d{3}_R[12]_001\.fastq\.gz

    Args:
        filename (str): The name of the file to check.
        sample_name (str): The sample name to use in the format check.

    Returns:
        bool: True if filename matches the 10x Genomics format, False otherwise.
    """

    pattern = re.compile(rf'^{re.escape(sample_name)}_S\d+_L\d{{3}}_R[12]_001\.fastq\.gz$')
    return pattern.match(filename)

def rename_fastqs(sample_name, fastq_dir):
    """
    Renames FASTQ files in the specified directory based on the sample name provided.

    Args:
        fastq_dir (str): Directory containing the FASTQ files.
        sample_name (str): Sample name to use in the new filename.
    Returns:
        string: Rename files
    """
    for filename in os.listdir(fastq_dir):
        if filename.endswith(".fastq.gz"):

            if not is_already_10x_format(sample_name, filename):

                print(f"Renaming {sample_name} {filename}")

                # Extract R1/R2 from _R1_ or _R2_ in the filename
                r_match = re.search(r'_R([12])_', filename)
                pair_number = r_match.group(1) if r_match else '1'

                # Extract lane number (e.g. L001, L002)
                l_match = re.search(r'_L(\d{3})_', filename)
                lane = l_match.group(1) if l_match else '001'

                # Extract sample number from _S<n>_L pattern to avoid matching sample name
                s_match = re.search(r'_S(\d+)_L', filename)
                sample_num = s_match.group(1) if s_match else '1'

                new_filename = f"{sample_name}_S{sample_num}_L{lane}_R{pair_number}_001.fastq.gz"

                old_path = os.path.join(fastq_dir, filename)
                new_path = os.path.join(fastq_dir, new_filename)
                os.rename(old_path, new_path)

            else:
                print(f"This is 10x ready {filename}")

if __name__ == "__main__":
    sample_name, fastq_dir = sys.argv[1], sys.argv[2]
    rename_fastqs(sample_name, fastq_dir)
