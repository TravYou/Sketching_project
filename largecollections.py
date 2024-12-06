import math
import mmh3
import time
import tracemalloc
from sourmash import MinHash, SourmashSignature
from sourmash.sbt import SBT, GraphFactory, Node, Leaf
from sourmash.sbtmh import SigLeaf
import subprocess
import os

# SBT Search Function Implementation (class SearchFunction created with assistance from Perplexity)
class SearchFunction:
    def __init__(self, threshold):
        self.threshold = threshold

    def __call__(self, node, query):
        return node.data.minhash.similarity(query.minhash) >= self.threshold

    def check_is_compatible(self, query):
        return True
    
    def score_fn(self, query_size, shared_size, subject_size, total_size):
        if query_size == 0 or subject_size == 0:
            return 0.0
        return shared_size / float(query_size)
    
    def passes(self, score):
        return score >= self.threshold
    
# Bloom Filter Implementation
class BloomFilter:
    def __init__(self, num_bits, num_hashes):
        self.num_bits = num_bits
        self.num_hashes = num_hashes
        self.bit_array = [0] * self.num_bits
        self.size = 0    

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

def brute_force_storage(bf_num_bits, capacity, dump_file_paths, query_kmer):
    total_storage = {}
    total_insertion_time = 0
    total_inserted_elements = 0
        
    tracemalloc.start()
    start_time = time.time()
    for file_path in dump_file_paths:
        new_bf, insertion_time = create_bloom_filter(bf_num_bits, capacity, file_path)
        total_insertion_time += insertion_time
        total_inserted_elements += new_bf.size
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        total_storage[base_name] = new_bf
        
    total_construction_time = time.time() - start_time
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\nQuerying for the K-mer...")
    tracemalloc.start()
    start_time = time.time()
    for file_name, filter in total_storage.items():
        if query_kmer in filter:
            print(f"{file_name} contains the query")
        else:
            print(f"{file_name} does not contain the query")

    print("Completed Querying for the K-mer...")
    total_query_time = time.time() - start_time
    current_query, peak_query = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Summary ---")
    print(f"Brute Force Storage Method:")
    print(f"  Inserted {total_inserted_elements} elements")
    print(f"  Insertion Time: {total_insertion_time:.2f} seconds")
    print(f"  Total Data Structure Construction Time: {total_construction_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")
    print(f"  Query Time: {total_query_time:.2f} seconds")
    print(f"  Querying Peak Memory Footprint: {peak_query / 1e6:.2f} MB")
        
def aggregate_bloom_filter(bf_num_bits, capacity, dump_file_paths, query_kmer):
    total_insertion_time = 0
    total_inserted_elements = 0
    
    tracemalloc.start()
    start_time = time.time()
    aggregated_filter = BloomFilter(bf_num_bits, int((bf_num_bits / capacity) * math.log(2)))

    for f_file in dump_file_paths:
        bf, insertion_time = create_bloom_filter(bf_num_bits, capacity, f_file)
        total_insertion_time += insertion_time
        total_inserted_elements += bf.size

        for bit_index in range(bf_num_bits):
            aggregated_filter.bit_array[bit_index] |= bf.bit_array[bit_index]
        
    aggregate_construction_time = time.time() - start_time
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\nQuerying for the K-mer...")
    tracemalloc.start()
    
    if query_kmer in aggregated_filter:
            print(f"Aggregate Filter contains the query")
    else:
        print(f"Aggregate Filter does not contain the query")
    
    print("Completed Querying for the K-mer...")
    total_query_time = time.time() - start_time
    current_query, peak_query = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Summary ---")
    print(f"Aggregate Bloom Filter Storage Method:")
    print(f"  Inserted {total_inserted_elements} elements")
    print(f"  Insertion Time: {total_insertion_time:.2f} seconds")
    print(f"  Aggregate Filter Construction Time: {aggregate_construction_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")
    print(f"  Query Time: {total_query_time:.2f} seconds")
    print(f"  Querying Peak Memory Footprint: {peak_query / 1e6:.2f} MB")

def create_signature(f_file, ksize, scaled, inserted_elements):
    minh = MinHash(n=0, ksize=ksize, scaled=scaled)
    start_time = time.time()
    
    with open(f_file, 'r') as file:
        for index, line in enumerate(file):
            if index % 2 == 1: 
                kmer = line.strip()     
                minh.add_sequence(kmer)
                inserted_elements += 1
                
    insertion_time = time.time() - start_time
    return SourmashSignature(minh), inserted_elements, insertion_time

def sequence_bloom_tree(fasta_files, query_kmer, table_size, n_tables):
    total_inserted_elements = 0
    total_insertion_time = 0
    ksize = len(query_kmer)
    scaled = 1000
    
    tracemalloc.start()
    start_time = time.time()
    factory = GraphFactory(ksize, table_size, n_tables)
    tree = SBT(factory, d=2)
    for f_file in fasta_files:
        base_name = os.path.splitext(os.path.basename(f_file))[0]
        sig, total_inserted_elements, insertion_time = create_signature(f_file, ksize, scaled, total_inserted_elements)
        total_insertion_time += insertion_time
        node = Node(factory, name=f"node_{base_name}")
        node.data.minhash = sig.minhash
        node.metadata['min_n_below'] = len(sig.minhash) 
        tree.add_node(node)
        
    total_construction_time = time.time() - start_time
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\nQuerying for the K-mer...")
    tracemalloc.start()
    start_time = time.time()
    query_mh = MinHash(n=0, ksize=ksize, scaled=scaled)
    query_mh.add_sequence(query_kmer)
    query_sig = SourmashSignature(query_mh)
    
    results = []
    for node in tree.find(search_fn=SearchFunction(1.0), query=query_sig):
        results.append(node.name)
    
    if results:
        print(f"SBT contains the query in the following files:")
        for name in results:
            print(f"- {name}")
    else:
        print("SBT does not contain the query")
        
    print("Completed Querying for the K-mer...")
    query_construction_time = time.time() - start_time
    current_query, peak_query = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print("\n--- Summary ---")
    print(f"SBT Storage Method:")
    print(f"  Inserted {total_inserted_elements} elements")
    print(f"  Insertion Time: {total_insertion_time:.2f} seconds")
    print(f"  SBT Construction Time: {total_construction_time:.2f} seconds")
    print(f"  Current Memory Usage: {current / 1e6:.2f} MB")
    print(f"  Peak Memory Footprint: {peak / 1e6:.2f} MB")
    print(f"  Query Time: {query_construction_time:.2f} seconds")
    print(f"  Querying Peak Memory Footprint: {peak_query / 1e6:.2f} MB")

def jellyfish_count_and_dump(input_file, dump_file, k_size, output_prefix="output"):
    count_cmd = f"jellyfish count -m {k_size} -s 100M -t 10 -C -o {output_prefix} {input_file}"
    subprocess.run(count_cmd, shell=True, check=True)

    dump_cmd = f"jellyfish dump {output_prefix} > {dump_file}"
    subprocess.run(dump_cmd, shell=True, check=True)

    if os.path.exists(dump_file):
        return os.path.abspath(dump_file)
    else:
        return None
    
def retrieve_kmer(file_path, kmer_num):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        if len(lines) >= 2:
            second_line_value = lines[1 + (kmer_num - 1) * 2].strip()
            return second_line_value
        
    return None
            
def main():
    # Total bits allocated for each filter
    TOTAL_BITS = 64 * 1024 * 1024 * 8  # 64 MB in bits
    table_size = 64 * 1024 * 1024
    
    print(f"Total bits allocated for each filter: {TOTAL_BITS}")

    # Target capacity: number of items
    capacity = 33_554_432

    # Bloom Filter Config
    bf_num_bits = TOTAL_BITS
    bf_bits_per_item = bf_num_bits / capacity
    bf_num_hashes = max(1, int(bf_bits_per_item * math.log(2)))
    n_tables = bf_num_hashes

    # Collections of reads files
    print("\nParsing Fasta Files...")
    fasta_file_paths = ["/Users/asmitha/CompGenomics/U00096_1.fasta", 
                        "/Users/asmitha/CompGenomics/U00096_2.fasta", 
                        "/Users/asmitha/CompGenomics/U00096_3.fasta",
                        "/Users/asmitha/CompGenomics/AE000509_1.fasta",
                        "/Users/asmitha/CompGenomics/AE000510_1.fasta"]
        
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
    
    # Retrieve Kmer to use for Querying
    print("\nRetrieving K-mer to use for Querying...")
    query_kmer = retrieve_kmer(dump_file_paths[0], 100)
    print(f"K-mer to query for is: {query_kmer}")
    
    # Brute Force storage of large collections of reads
    print("\nRunning Brute Force Storage Method...")
    brute_force_storage(bf_num_bits, capacity, dump_file_paths, query_kmer)
    
    # Aggregate Bloom Filter
    print("\nRunning Combined Bloom Filter Storage Method...")
    aggregate_bloom_filter(bf_num_bits, capacity, dump_file_paths, query_kmer)
    
    # Sequence Bloom Tree
    print("\nRunning Sequence Bloom Tree Storage Method...")
    sequence_bloom_tree(dump_file_paths, query_kmer, table_size, n_tables)

if __name__ == "__main__":
    main()
