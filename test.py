import sys
import os
import subprocess
import networkx as nx
from Bio import SeqIO

current_script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_script_path, 'Query'))
from read_and_query import query_graph

def generate_mock_data():
    subprocess.run(["python3", os.path.join(current_script_path, "Datasets", "generate_mock_genomic_data.py")])

def build_wheeler_graph(fasta_files):
    args = ["python3", os.path.join(current_script_path, "Build_WG", "build_debruijn.py"), "-o", "wheeler_graph.dot"] + fasta_files
    subprocess.run(args)

def query_reads(graph_file, fastq_file):
    G = nx.nx_pydot.read_dot(graph_file)
    predictions = {}
    with open(fastq_file, 'r') as f:
        for record in SeqIO.parse(f, "fastq"):
            read = str(record.seq)
            organism_ids = query_graph(G, read)
            predictions[record.id] = organism_ids  # Store predictions by read ID
    return predictions

def compare_with_ground_truth(predictions, fastq_file):
    correct_predictions = 0
    total_predictions = 0

    with open(fastq_file, 'r') as f:
        for record in SeqIO.parse(f, "fastq"):
            # Extract ground truth from the header
            header_parts = record.id.split('|')
            source_label = header_parts[0]
            true_organism_id = {"host": 0, "contaminant1": 1, "contaminant2": 2}.get(source_label, -1)

            # Compare with prediction
            predicted_organism_ids = predictions.get(record.id, set())
            if true_organism_id in predicted_organism_ids:
                correct_predictions += 1
            total_predictions += 1

    accuracy = correct_predictions / total_predictions
    print(f"Accuracy: {accuracy:.2f}")

if __name__ == "__main__":
    generate_mock_data()
    fasta_files = ['host.fasta', 'contaminant1.fasta', 'contaminant2.fasta']
    build_wheeler_graph(fasta_files)
    fastq_file = 'mock_data.fastq'
    predictions = query_reads("wheeler_graph.dot", fastq_file)
    compare_with_ground_truth(predictions, fastq_file)

