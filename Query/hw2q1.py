import sys

# Ensure valid fatsq file
def is_valid_input_file(input):
    # Every valid FASTQ record should have exactly 4 lines
    if len(input) % 4 != 0:
        return False, "Error: The FASTQ file does not have a multiple of 4 lines."

    for i in range(0, len(input), 4):
        # The first line of each record should start with the "@" symbol
        if not input[i].startswith("@"):
            return False, f"Error: Line {i + 1} does not start with '@'."
        
        # The third line of each record should start with the "+" symbol
        if not input[i + 2].startswith("+"):
            return False, f"Error: Line {i + 3} does not start with '+'."
        
        # The sequence and quality strings should have the same length
        if len(input[i + 1].strip()) != len(input[i + 3].strip()):
            return False, f"Error: Length mismatch between sequence and quality on lines {i + 2} and {i + 4}."
    
    return True, ""

# Function to read input file
def read_input_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Make sure we have a valid input file
    is_valid, error_msg = is_valid_input_file(lines)
    if not is_valid:
        print(error_msg)
        sys.exit(1)

    sequences = lines[1::4]  # Starting on line 1 every four lines
    quality_chars = lines[3::4]  # Starting on line 3 every fout lines
    quality_scores = [[ord(char) - 33 for char in line.strip()] for line in quality_chars]

    return sequences, quality_scores

# Function to write to output file
def output_to_desired_file(filename, output):
    with open(filename, 'w') as file:
        file.write(" ".join(map(str, output)))

# Total quality
def get_total_quality(quality_scores):
    total_quality = [sum(scores) for scores in quality_scores]
    return total_quality

# The read that has the lowest total quality score (1=first read, 2=second read, etc)
def get_lowest_total_quality(total_quality):
    return total_quality.index(min(total_quality)) + 1


# The read that has the highest total quality score
def get_highest_total_quality(total_quality):
    return total_quality.index(max(total_quality)) + 1

# The total number of base quality values (denoted Q in lecture) less than 10
def get_total_Q_less_than_10(quality_scores):
    return sum([1 for read in quality_scores for q in read if q < 10])

# The total number of base quality values ≥ 30
def get_total_Q_greater_than_30(quality_scores):
    return sum([1 for read in quality_scores for q in read if q >= 30])

# Total number of characters observed in the read sequences of the file other than A, C, G, or T
def total_unknown_bases(sequences):
    return sum([1 for read in sequences for base in read.strip() if base not in ["A", "C", "G", "T"]])


def main():
    if len(sys.argv) != 3:
        print("Incorrect argument usage. Make sure it's in the following form: python hw2q1.py <input_filename> <output_filename>")
        return
    
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    sequences, quality_scores = read_input_file(input_filename)
    total_quality = get_total_quality(quality_scores)
    
    output = [
        get_lowest_total_quality(total_quality),
        get_highest_total_quality(total_quality),
        get_total_Q_less_than_10(quality_scores),
        get_total_Q_greater_than_30(quality_scores),
        total_unknown_bases(sequences)
    ]

    output_to_desired_file(output_filename, output)


if __name__ == "__main__":
    main()