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

# Bloom Filter Implementation
class BloomFilter:
    def __init__(self, num_bits, num_hashes):
        self.num_bits = num_bits
        self.num_hashes = num_hashes
        self.bit_array = [0] * self.num_bits
        self.size = 0  # Number of inserted elements

    def add(self, element):
        for i in range(self.num_hashes):
            hash_val = mmhmm.hash(element, i, signed=False) % self.num_bits
            self.bit_array[hash_val] = 1
        self.size += 1

    def __contains__(self, element):
        for i in range(self.num_hashes):
            hash_val = mmhmm.hash(element, i, signed=False) % self.num_bits
            if self.bit_array[hash_val] == 0:
                return False
        return True

# Function to generate random strings
def generate_random_strings(num_strings, string_length=10):
    random_strings = set()
    while len(random_strings) < num_strings:
        s = ''.join(random.choices(string.ascii_letters + string.digits, k=string_length))
        random_strings.add(s)
    return list(random_strings)

# Function to measure false positive rate
def measure_false_positive_rate(filter_obj, filter_name, test_data, num_test_items):
    false_positives = 0
    for item in test_data:
        if item in filter_obj:
            false_positives += 1
    false_positive_rate = false_positives / num_test_items
    print(f"{filter_name} False Positive Rate: {false_positive_rate * 100:.4f}%")
    return false_positive_rate

def create_bloom_filter(bf_num_bits, capacity, f_file):
    bf_bits_per_item = bf_num_bits / capacity
    bf_num_hashes = max(1, int(bf_bits_per_item * math.log(2)))
    bf = BloomFilter(num_bits=bf_num_bits, num_hashes=bf_num_hashes)
    
    print(f"\n--- Bloom Filter Configuration ---")
    print(f"Number of bits: {bf_num_bits}")
    print(f"Number of hashes: {bf_num_hashes}")
    print(f"Capacity: {capacity} items")
    print(f"Bits per item: {bf_bits_per_item:.2f} bits")
    
    # Generate data for insertion
    k_len = 10
    print(f"\nGenerating k-mers with k={k_len} for insertion...")
    data = retrieve_kmers(k_len, f_file)
    print(len(data))

    # Start memory tracking
    tracemalloc.start()
    
    # Insert elements into Bloom Filter
    print("\nInserting into Bloom Filter...")
    start_time = time.time()
    for key, value in data.items():
        bf.add(key)
    bf_insertion_time = time.time() - start_time
    print(f"Bloom Filter inserted {bf.size} elements in {bf_insertion_time:.2f} seconds.")
    
    # Stop memory tracking
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    # Prepare test data for false positive rate measurement
    num_test_items = 100000
    test_data = generate_random_strings(num_test_items)

    print("\nMeasuring False Positive Rates...")
    bf_fpr = measure_false_positive_rate(bf, "Bloom Filter", test_data, num_test_items)
    
    print("\n--- Summary ---")
    print(f"Bloom Filter:")
    print(f"  Inserted {bf.size} elements")
    print(f"  Bits per item: {bf_bits_per_item:.2f}")
    print(f"  False Positive Rate: {bf_fpr * 100:.4f}%")
    print(f"  Construction Time: {bf_insertion_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")
    
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


def brute_force_storage(bf_num_bits, capacity, fasta_files):
    total_storage = {}
    total_insertion_time = 0
    total_inserted_elements = 0
    
    tracemalloc.start()
    
    for f_file, bases in fasta_files.items():
        new_bf, insertion_time = create_bloom_filter(bf_num_bits, capacity, bases)
        total_insertion_time += insertion_time
        total_inserted_elements += new_bf.size
        total_storage[i] = new_bf
        
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
        
    print("\n--- Summary ---")
    print(f"Brute Force:")
    print(f"  Inserted {total_inserted_elements} elements")
    print(f"  Construction Time: {total_insertion_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")
    
def aggregate_bloom_filter(bf_num_bits, capacity, fasta_files):
    aggregated_filter = BloomFilter(f_num_bits, int((f_num_bits / capacity) * math.log(2)))
    total_insertion_time = 0

    tracemalloc.start(bf_num_bits, capacity, fasta_file_paths)

    for f_file in fasta_files:
        bf, insertion_time = create_bloom_filter(bf_num_bits, capacity, f_file)
        total_insertion_time += insertion_time

        # Aggregate current file's Bloom Filter into the combined filter
        for bit_index in range(f_num_bits):
            aggregated_filter.bit_array[bit_index] |= bf.bit_array[bit_index]

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Aggregative Summary ---")
    print(f"Total Construction Time: {total_insertion_time:.2f} seconds")
    print(f"Combined Bloom Filter Memory Usage: {current / 1e6:.2f} MB")
    print(f"Peak Memory Footprint: {peak / 1e6:.2f} MB")
    
def sequence_bloom_tree(bf_num_bits, capacity, fasta_files):
    
    tracemalloc.start()
    
    for f_file in fasta_files:
        bloom_filter = create_bloom_filter(bf_num_bits, capacity, f_file)
        sbt_filename = "bloom_filter_sbt.sbt.zip"
        insert_bloom_filter_into_sbt(sbt_filename, bloom_filter, name="bloom_filter_sbt")
        
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Sequence Bloom Filter Summary ---")
    print(f"Total Construction Time: {total_insertion_time:.2f} seconds")
    print(f"Combined Bloom Filter Memory Usage: {current / 1e6:.2f} MB")
    print(f"Peak Memory Footprint: {peak / 1e6:.2f} MB")
    
def read_fasta_files(file_path):
    print(file_path)
    with open(file_path, 'r') as in_file:
        final_line = ""
        for line in in_file:
            if line[0] == '>':
                continue
            else:
                line=line.strip()
                for c in line:
                    if c in 'ACGT':
                        final_line = final_line + c
    print(len(final_line))
    return final_line

def retrieve_kmers(k_len, fasta_file):
    print(len(fasta_file))
    index = {}
    start = 0
    end = int(k_len)
    while end < len(fasta_file)+1:
        sub_string = str(fasta_file[start:end])
        if sub_string not in index:
                index[sub_string] = [start]
        else:
            index[sub_string].append(start)

        end += 1
    return index
    
def main():
    # Total bits allocated for each filter
    TOTAL_BITS = 64 * 1024 * 1024 * 8  # 64 MB in bits
    print(f"Total bits allocated for each filter: {TOTAL_BITS}")

    # Target capacity: number of items
    capacity = 33_554_432

    # Bloom Filter Config
    bf_num_bits = TOTAL_BITS
    
    # Collections of reads files
    print("\nParsing Fasta Files...")
    fasta_file_paths = ["U00096_1.fasta", "U00096_2.fasta", "U00096_3.fasta"]
    fasta_files = {}
    for f_file in fasta_file_paths:
        base_name = f_file.split('.')[0]
        fasta_files[base_name] = read_fasta_files(f_file)

    # Brute Force storage of large collections of reads
    print("\nRunning Brute Force Storage Method...")
    brute_force_storage(bf_num_bits, capacity, fasta_files)
    
    # Aggregate Bloom Filter
    print("\nRunning Combined Bloom Filter Storage Method...")
    aggregate_bloom_filter(bf_num_bits, capacity, fasta_files)
    
    # Sequence Bloom Tree
    print("\nRunning Sequence Bloom Tree Storage Method...")
    sequence_bloom_tree(bf_num_bits, capacity, fasta_files)
    
    

if __name__ == "__main__":
    main()
