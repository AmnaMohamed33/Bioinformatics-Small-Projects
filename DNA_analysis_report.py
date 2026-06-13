#------------ PROJECT's STEPS ----------
#  download a genome from NCBI
# calculate the GC content
# enumerate the nucleotides
# convert RNA to protein
# print a report
# -----------------------------------------------------

import Bio                        # Import the main Biopython library
from Bio.Seq import Seq           # Import Seq object for sequence manipulation
from Bio import SeqIO             # Import SeqIO for reading and writing sequence files

max_gc=0                          # Initialize maximum GC content as 0

# Loop through all sequences in the FASTA file
for record in SeqIO.parse("/home/useramna/ROSALIND-Examples.py/sequence.fasta", "fasta"):

    # Calculate GC content by counting G and C bases divided by total length multiplied by 100
    gc=(record.seq.count('G') + record.seq.count('C')) / (len(record.seq)) *100
    if gc > max_gc:                                # Check if current GC is higher than the maximum found so far
        max_gc = gc                                # Update maximum GC content

# Print the maximum GC content rounded to 3 decimal places
print(f"the GC content percentage of this sequence is: {max_gc:.3f} %")         
seq=record.seq                             # Store the last sequence from the loop

# Count and print the total number of nucleotides in the sequence
print(f"the number of total nucleotides of this sequence is: {seq.count("A")+ seq.count("T")+ seq.count("G") + seq.count("C")} Nucleotides")

rna=seq.transcribe()                             # Transcribe DNA sequence to RNA by replacing T with U
print(f"the proteins translated from this sequence are: {rna.translate()}")               # Translate RNA to protein and print the result


# Write the analysis results to a text report file
with open("DNA_analysis_report.txt", "w")as f:
    f.write("DNA Analysis Report\n")                    # Write report title
    f.write(f"GC content: {max_gc:.3f}% \n")            # Write GC content
    f.write(f"Total Nucleotides :{seq.count("A")+ seq.count("T")+ seq.count("G") + seq.count("C")} Nucleotides\n")        # Write total nucleotide count
    f.write(f"RNA: {rna}\n")                            # Write RNA sequence
    f.write(f"Protein:{rna.translate()}\n")             # Write protein sequence




