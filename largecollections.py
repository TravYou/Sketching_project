import math
import mmh3
import random
import string
import sys
import gc
import time
import tracemalloc
from sourmash import MinHash, SourmashSignature
from sourmash.sbt import SBT, Node
from sourmash.sbtmh import SigLeaf
import subprocess
import os



# Bloom Filter Implementation
class BloomFilter:
    def __init__(self, num_bits, num_hashes):
        self.num_bits = num_bits
        self.num_hashes = num_hashes
        self.bit_array = [0] * self.num_bits
        self.size = 0  # Number of inserted elements

    def add(self, element):
        for i in range(self.num_hashes):
            hash_val = mmh3.hash(element, i, signed=False) % self.num_bits
            self.bit_array[hash_val] = 1
        self.size += 1

    def __contains__(self, element):
        for i in range(self.num_hashes):
            hash_val = mmh3.hash(element, i, signed=False) % self.num_bits
            if self.bit_array[hash_val] == 0:
                return False
        return True

def create_bloom_filter(bf_num_bits, capacity, f_file):
    bf_bits_per_item = bf_num_bits / capacity
    bf_num_hashes = max(1, int(bf_bits_per_item * math.log(2)))
    bf = BloomFilter(num_bits=bf_num_bits, num_hashes=bf_num_hashes)
    
    # Insert elements into Bloom Filter
    print("\nInserting into Bloom Filter...")
    start_time = time.time()
    
    with open(f_file, 'r') as file:
        for index, line in enumerate(file):
            if index % 2 == 1: 
                kmer = line.strip()
                bf.add(kmer)
    print("Completed Inserting Elements...")
    bf_insertion_time = time.time() - start_time
    print(f"Bloom Filter inserted {bf.size} elements in {bf_insertion_time:.2f} seconds.")
    
    return bf, bf_insertion_time
    
def insert_bloom_filter_into_sbt(sbt_file, bloom_filter, name):
    try:
        sbt = SBT.load(sbt_file)
    except FileNotFoundError:
        created_nodes = GraphFactory(ksize=bloom_filter.ksize(), starting_size=bloom_filter.size(), n_tables=bloom_filter.n_tables())
        sbt = SBT(created_nodes, d=2)

    leaf = SigLeaf(name, bloom_filter)
    sbt.add_node(leaf)

    sbt.save(sbt_file)

def brute_force_storage(bf_num_bits, capacity, dump_file_paths):
    total_storage = {}
    total_insertion_time = 0
    total_inserted_elements = 0
        
    for file_path in dump_file_paths:
        new_bf, insertion_time = create_bloom_filter(bf_num_bits, capacity, file_path)
        total_insertion_time += insertion_time
        total_inserted_elements += new_bf.size
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        total_storage[base_name] = new_bf
    
    return total_inserted_elements, total_insertion_time
    
def aggregate_bloom_filter(bf_num_bits, capacity, dump_file_paths):
    aggregated_filter = BloomFilter(bf_num_bits, int((bf_num_bits / capacity) * math.log(2)))
    total_insertion_time = 0
    total_inserted_elements = 0

    for f_file in dump_file_paths:
        bf, insertion_time = create_bloom_filter(bf_num_bits, capacity, f_file)
        total_insertion_time += insertion_time
        total_inserted_elements += bf.size

        # Aggregate current file's Bloom Filter into the combined filter
        for bit_index in range(bf_num_bits):
            aggregated_filter.bit_array[bit_index] |= bf.bit_array[bit_index]
    
    return total_inserted_elements, total_insertion_time
    
def sequence_bloom_tree(bf_num_bits, capacity, fasta_files):
        
    for f_file in fasta_files:
        bloom_filter = create_bloom_filter(bf_num_bits, capacity, f_file)
        sbt_filename = "bloom_filter_sbt.sbt.zip"
        insert_bloom_filter_into_sbt(sbt_filename, bloom_filter, name="bloom_filter_sbt")
        
def jellyfish_count_and_dump(input_file, dump_file, k_size, output_prefix="output"):
    # Count k-mers
    count_cmd = f"jellyfish count -m {k_size} -s 100M -t 10 -C -o {output_prefix} {input_file}"
    subprocess.run(count_cmd, shell=True, check=True)

    # Dump k-mers to text file
    dump_cmd = f"jellyfish dump {output_prefix} > {dump_file}"
    subprocess.run(dump_cmd, shell=True, check=True)

    # Check if the dump file was created and return its path
    if os.path.exists(dump_file):
        return os.path.abspath(dump_file)
    else:
        return None
    
def main():
    # Total bits allocated for each filter
    TOTAL_BITS = 64 * 1024 * 1024 * 8  # 64 MB in bits
    print(f"Total bits allocated for each filter: {TOTAL_BITS}")
    
    # Log file for writing analysis to
    brute_force_log_file = "brute_force_log_file.log"
    aggregate_bloom_filter_log_file = "aggregate_bloom_filter.log"
    sequence_bloom_tree_log_file = "sequence_bloom_tree_log_file.log"
    

    # Target capacity: number of items
    capacity = 33_554_432

    # Bloom Filter Config
    bf_num_bits = TOTAL_BITS
    
    # Collections of reads files
    print("\nParsing Fasta Files...")
    fasta_file_paths = ["/Users/asmitha/CompGenomics/U00096_1.fasta", 
                        "/Users/asmitha/CompGenomics/U00096_2.fasta", 
                        "/Users/asmitha/CompGenomics/U00096_3.fasta"]
        
    # Creating list of k-mers
    print("\nCreating Lists of K-mers...")
    k_len = 100
    dump_file_paths = []
    total_construction_time = []
    
    tracemalloc.start()
    start_time = time.time()
    for file_path in fasta_file_paths:
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        dump_file = f"{base_name}_kmers.txt"
    
        print(f"\nGenerating k-mers for {base_name} with k={k_len}...")
        dump_file_path = jellyfish_count_and_dump(file_path, dump_file, k_len)
        
        dump_file_paths.append(dump_file_path)

    total_construction_time = time.time() - start_time
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Summary ---")
    print(f"Creating K-Mers:")
    print(f"  Construction Time: {total_construction_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")
    
    # Brute Force storage of large collections of reads
    print("\nRunning Brute Force Storage Method...")
    
    tracemalloc.start()
    total_inserted_elements, total_insertion_time = brute_force_storage(bf_num_bits, capacity, dump_file_paths)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Summary ---")
    print(f"Brute Force Storage Method:")
    print(f"  Inserted {total_inserted_elements} elements")
    print(f"  Construction Time: {total_insertion_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")
    
    # Aggregate Bloom Filter
    print("\nRunning Combined Bloom Filter Storage Method...")
    
    tracemalloc.start()
    total_inserted_elements, total_insertion_time = aggregate_bloom_filter(bf_num_bits, capacity, dump_file_paths)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Summary ---")
    print(f"Aggregated Bloom Filter Storage Method:")
    print(f"  Inserted {total_inserted_elements} elements")
    print(f"  Construction Time: {total_insertion_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")
    
    # Sequence Bloom Tree
    print("\nRunning Sequence Bloom Tree Storage Method...")
    
    tracemalloc.start()
    sequence_bloom_tree(bf_num_bits, capacity, dump_file_paths)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Summary ---")
    print(f"Sequence Bloom Filter Storage Method:")
    print(f"  Inserted {total_inserted_elements} elements")
    print(f"  Construction Time: {total_insertion_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")  
    
    

if __name__ == "__main__":
    main()

#conda create -n bioenv python=3.12
#conda activate bioenv
#conda install -c conda-forge sourmash-minimal
#conda install mmh3
#conda install -c conda-forge biopython
#python3 largecollections.py 