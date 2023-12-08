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
    recognized_reads_count = 0
    with open(fastq_file, 'r') as f:
        for record in SeqIO.parse(f, "fastq"):
            read = str(record.seq)
            organism_ids = query_graph(G, read)
            predictions[record.id] = organism_ids
            if organism_ids:
                recognized_reads_count += 1
    return predictions, recognized_reads_count

def compare_with_ground_truth(predictions, fastq_file):
    true_positives = {0: 0, 1: 0, 2: 0}
    false_positives = {0: 0, 1: 0, 2: 0}
    false_negatives = {0: 0, 1: 0, 2: 0}
    true_negatives = {0: 0, 1: 0, 2: 0}
    ground_truth_counts = {0: 0, 1: 0, 2: 0}
    predicted_counts = {0: 0, 1: 0, 2: 0}
    total_reads = 0
    ambiguous = 0

    with open(fastq_file, 'r') as f:
        for record in SeqIO.parse(f, "fastq"):
            total_reads += 1
            header_parts = record.id.split('|')
            source_label = header_parts[0]
            true_organism_id = {"host": 0, "contaminant1": 1, "contaminant2": 2}.get(source_label, -1)
            ground_truth_counts[true_organism_id] += 1

            predicted_organism_ids = predictions.get(record.id, set())

            if len(predicted_organism_ids) > 1:
                ambiguous += 1
                continue

            predicted_organism_id = next(iter(predicted_organism_ids), -1)
            predicted_counts[predicted_organism_id] += 1

            if true_organism_id == predicted_organism_id:
                true_positives[true_organism_id] += 1
            else:
                false_negatives[true_organism_id] += 1
                if predicted_organism_id != -1:
                    false_positives[predicted_organism_id] += 1

    for org_id in [0, 1, 2]:
        true_negatives[org_id] = total_reads - (true_positives[org_id] + false_positives[org_id] + false_negatives[org_id])

    metrics = {}
    for org_id in [0, 1, 2]:
        precision = true_positives[org_id] / (true_positives[org_id] + false_positives[org_id]) if true_positives[org_id] + false_positives[org_id] > 0 else 0
        recall = true_positives[org_id] / (true_positives[org_id] + false_negatives[org_id]) if true_positives[org_id] + false_negatives[org_id] > 0 else 0
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        accuracy = (true_positives[org_id] + true_negatives[org_id]) / total_reads if total_reads > 0 else 0

        metrics[org_id] = {"precision": precision, "recall": recall, "f1_score": f1_score, "accuracy": accuracy}

    for org_id in [0, 1, 2]:
        category = {0: "Host", 1: "Contaminant1", 2: "Contaminant2"}[org_id]
        print(f"Total Ground Truth for {category}: {ground_truth_counts[org_id]}")
        print(f"Total Predicted for {category}: {predicted_counts[org_id]}")
        print(f"Metrics for {category} - Precision: {metrics[org_id]['precision']:.2f}, Recall: {metrics[org_id]['recall']:.2f}, F1 Score: {metrics[org_id]['f1_score']:.2f}, Accuracy: {metrics[org_id]['accuracy']:.2f}")
    print(f"Ambiguous Classifications: {ambiguous}")

if __name__ == "__main__":
    generate_mock_data()
    fasta_files = ['host.fasta', 'contaminant1.fasta', 'contaminant2.fasta']
    build_wheeler_graph(fasta_files)
    fastq_file = 'mock_data.fastq'
    predictions, recognized_reads_count = query_reads("wheeler_graph.dot", fastq_file)
    print(f"Number of reads recognized by the graph: {recognized_reads_count}")
    compare_with_ground_truth(predictions, fastq_file)