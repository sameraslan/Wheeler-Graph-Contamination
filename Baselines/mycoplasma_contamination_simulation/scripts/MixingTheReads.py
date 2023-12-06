import random

def mix_reads(source_file1, source_file2, target_file, contamination_rate):
    with open(source_file1, 'r') as f1, open(source_file2, 'r') as f2, open(target_file, 'w') as out:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        
        mycoplasma_contaminated_reads = int(len(lines1) * contamination_rate / 4) * 4
        
        for i in range(0, len(lines1), 4):
            if mycoplasma_contaminated_reads > 0 and random.random() < contamination_rate:
                if len(lines2) > i:  # Check if there's a corresponding read in Mycoplasma reads
                    out.writelines(lines2[i:i+4])
                    mycoplasma_contaminated_reads -= 4
                else:
                    out.writelines(lines1[i:i+4])  # Write human read if no Mycoplasma read is available
            else:
                out.writelines(lines1[i:i+4])

# Example usage
mix_reads('human_single_end_reads.fq', 'mycoplasma_reads.fq', 'mixed_reads.fq', 0.05)
