#!/usr/bin/env python3

####!/opt/conda/bin/python

import os
import sys
import re

def is_already_10x_format(sample_name, filename):
    """
    Determine if the filename matches the expected 10x Genomics format.
    {sample_name}_S<n>_L<lane>_R<1|2>_001.fastq.gz

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

                # Read number from _R1_/_R2_, and also _R1./_R2. and _1./_2. endings
                # (X_R1.fastq.gz, fetchngs' SRR_1.fastq.gz). The original matched only
                # _R1_ and defaulted to 1, so a pair named X_R1.fastq.gz / X_R2.fastq.gz
                # BOTH became ..._R1_001 and R2 overwrote R1 -- half the reads gone, exit 0.
                # The LAST read token wins: the read number sits at the end of the name,
                # and a sample name may itself contain _R1_ (pt_R1_tcr_S1_L001_R2_001).
                # Never guess: an unreadable read number is an error.
                r_tokens = re.findall(r'_R([12])(?=[_.])', filename)
                if r_tokens:
                    pair_number = r_tokens[-1]
                else:
                    end_match = re.search(r'[_.]([12])\.f(?:ast)?q\.gz$', filename)
                    if not end_match:
                        sys.exit(f"ERROR: cannot tell R1 from R2 in '{filename}'. Name it "
                                 f"..._R1_001.fastq.gz / ..._R2_001.fastq.gz (or _R1/_R2).")
                    pair_number = end_match.group(1)

                # Extract lane number (e.g. L001, L002)
                l_match = re.search(r'_L(\d{3})_', filename)
                lane = l_match.group(1) if l_match else '001'

                # Extract sample number from _S<n>_L pattern to avoid matching sample name
                s_match = re.search(r'_S(\d+)_L', filename)
                sample_num = s_match.group(1) if s_match else '1'

                new_filename = f"{sample_name}_S{sample_num}_L{lane}_R{pair_number}_001.fastq.gz"

                old_path = os.path.join(fastq_dir, filename)
                new_path = os.path.join(fastq_dir, new_filename)
                # Two inputs mapping to one name (e.g. lanes from different flowcells both
                # L001) would silently drop a file. Refuse instead.
                if os.path.lexists(new_path):
                    sys.exit(f"ERROR: renaming '{filename}' -> '{new_filename}' would overwrite "
                             f"an existing file. Give the FASTQs distinct lane numbers.")
                os.rename(old_path, new_path)

            else:
                print(f"This is 10x ready {filename}")

if __name__ == "__main__":
    sample_name, fastq_dir = sys.argv[1], sys.argv[2]
    rename_fastqs(sample_name, fastq_dir)
