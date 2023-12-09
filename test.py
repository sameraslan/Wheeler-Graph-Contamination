import sys
import os
import subprocess
import networkx as nx
from Bio import SeqIO
from tabulate import tabulate

# Determine the path of the current script and add the 'Query' directory to the system path.
current_script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_script_path, 'Query'))

# Import the query_graph function from the read_and_query module.
from read_and_query import query_graph

# Function to generate mock genomic data by running an external Python script.
def generate_mock_data():
    subprocess.run(["python3", os.path.join(current_script_path, "Datasets", "generate_mock_genomic_data.py")])

# Function to build a Wheeler graph from a list of FASTA files.
def build_wheeler_graph(fasta_files):
    # Combine the command and arguments as a list, and execute it.
    args = ["python3", os.path.join(current_script_path, "Build_WG", "build_debruijn.py"), "-o", "wheeler_graph.dot"] + fasta_files
    subprocess.run(args)

# Function to query reads from a FASTQ file against a Wheeler graph.
def query_reads(graph_file, fastq_file):
    # Load the graph from a DOT file.
    G = nx.nx_pydot.read_dot(graph_file)
    predictions = {}
    recognized_reads_count = 0

    # Open the FASTQ file and parse each record.
    with open(fastq_file, 'r') as f:
        for record in SeqIO.parse(f, "fastq"):
            # Extract the sequence from the record and query the graph.
            read = str(record.seq)
            organism_ids = query_graph(G, read)
            predictions[record.id] = organism_ids

            # If the organism IDs are found, increment the recognized reads count.
            if organism_ids:
                recognized_reads_count += 1
    return predictions, recognized_reads_count

# Function to compare the predictions with ground truth data.
def compare_with_ground_truth(predictions, fastq_file):
    # Initialize counters for various metrics.
    true_positives = {0: 0, 1: 0, 2: 0}
    false_positives = {0: 0, 1: 0, 2: 0}
    false_negatives = {0: 0, 1: 0, 2: 0}
    true_negatives = {0: 0, 1: 0, 2: 0}
    ground_truth_counts = {0: 0, 1: 0, 2: 0}
    predicted_counts = {0: 0, 1: 0, 2: 0}
    total_reads = 0
    ambiguous = 0

    # Open the FASTQ file and parse each record for comparison.
    with open(fastq_file, 'r') as f:
        for record in SeqIO.parse(f, "fastq"):
            total_reads += 1
            header_parts = record.id.split('|')
            source_label = header_parts[0]
            true_organism_id = {"host": 0, "contaminant1": 1, "contaminant2": 2}.get(source_label, -1)
            ground_truth_counts[true_organism_id] += 1

            predicted_organism_ids = predictions.get(record.id, set())

            # Handle ambiguous classifications.
            if len(predicted_organism_ids) > 1:
                ambiguous += 1
                continue

            predicted_organism_id = next(iter(predicted_organism_ids), -1)
            predicted_counts[predicted_organism_id] += 1

            # Calculate true positives, false positives, and false negatives.
            if true_organism_id == predicted_organism_id:
                true_positives[true_organism_id] += 1
            else:
                false_negatives[true_organism_id] += 1
                if predicted_organism_id != -1:
                    false_positives[predicted_organism_id] += 1

    # Calculate true negatives and other metrics.
    for org_id in [0, 1, 2]:
        true_negatives[org_id] = total_reads - (true_positives[org_id] + false_positives[org_id] + false_negatives[org_id])

    metrics = {}
    for org_id in [0, 1, 2]:
        # Calculate precision, recall, F1 score, and accuracy.
        precision = true_positives[org_id] / (true_positives[org_id] + false_positives[org_id]) if true_positives[org_id] + false_positives[org_id] > 0 else 0
        recall = true_positives[org_id] / (true_positives[org_id] + false_negatives[org_id]) if true_positives[org_id] + false_negatives[org_id] > 0 else 0
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        accuracy = (true_positives[org_id] + true_negatives[org_id]) / total_reads if total_reads > 0 else 0

        metrics[org_id] = {"precision": precision, "recall": recall, "f1_score": f1_score, "accuracy": accuracy}

    # Print out the comparison results and metrics.
    for org_id in [0, 1, 2]:
        category = {0: "Host", 1: "Contaminant1", 2: "Contaminant2"}[org_id]
        print(f"Total Ground Truth for {category}: {ground_truth_counts[org_id]}")
        print(f"Total Predicted for {category}: {predicted_counts[org_id]}")
        print(f"Metrics for {category} - Precision: {metrics[org_id]['precision']:.2f}, Recall: {metrics[org_id]['recall']:.2f}, F1 Score: {metrics[org_id]['f1_score']:.2f}, Accuracy: {metrics[org_id]['accuracy']:.2f}")
    print(f"Ambiguous Classifications: {ambiguous}")

    # Function to display a summary table of the results.
    display_results_table(ground_truth_counts, predicted_counts, metrics, ambiguous)

# Function to display a formatted table of the results.
def display_results_table(ground_truth_counts, predicted_counts, metrics, ambiguous):
    table_data = []
    for org_id, category in {0: "Host", 1: "Contaminant1", 2: "Contaminant2"}.items():
        row = [category, 
               ground_truth_counts[org_id], 
               predicted_counts[org_id], 
               f"{metrics[org_id]['precision']:.2f}", 
               f"{metrics[org_id]['recall']:.2f}", 
               f"{metrics[org_id]['f1_score']:.2f}", 
               f"{metrics[org_id]['accuracy']:.2f}"]
        table_data.append(row)

    table_data.append(["Ambiguous", "0", ambiguous, "-", "-", "-", "-"])

    headers = ["Category", "Total Ground Truth", "Total Predicted", "Precision", "Recall", "F1 Score", "Accuracy"]
    table = tabulate(table_data, headers=headers, tablefmt="pretty")
    print(table)

# Main function to run the entire pipeline.
if __name__ == "__main__":
    generate_mock_data()
    fasta_files = ['host.fasta', 'contaminant1.fasta', 'contaminant2.fasta']
    build_wheeler_graph(fasta_files)
    fastq_file = 'mock_data.fastq'
    predictions, recognized_reads_count = query_reads("wheeler_graph.dot", fastq_file)
    print(f"Number of reads recognized by the graph: {recognized_reads_count}")
    compare_with_ground_truth(predictions, fastq_file)