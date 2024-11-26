import math
import mmh3
import random
import string
import sys
import gc
import time

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

# Cuckoo Filter Implementation
class CuckooFilter:
    def __init__(self, num_buckets, bucket_size, fingerprint_size, max_kicks=500):
        self.num_buckets = num_buckets
        self.bucket_size = bucket_size
        self.fingerprint_size = fingerprint_size
        self.max_kicks = max_kicks
        self.buckets = [[] for _ in range(self.num_buckets)]
        self.size = 0  # Number of inserted elements

    def _hash(self, item):
        return mmh3.hash(item, signed=False) % self.num_buckets

    def _fingerprint(self, item):
        fp = mmh3.hash(item, seed=42, signed=False) & ((1 << self.fingerprint_size) - 1)
        if fp == 0:
            fp += 1
        return fp

    def _alt_index(self, index, fp):
        fp_hash = mmh3.hash(str(fp), seed=43, signed=False)
        return (index ^ fp_hash) % self.num_buckets

    def add(self, item):
        fp = self._fingerprint(item)
        i1 = self._hash(item)
        i2 = self._alt_index(i1, fp)

        if len(self.buckets[i1]) < self.bucket_size:
            self.buckets[i1].append(fp)
            self.size += 1
            return True

        if len(self.buckets[i2]) < self.bucket_size:
            self.buckets[i2].append(fp)
            self.size += 1
            return True

        i = random.choice([i1, i2])
        for _ in range(self.max_kicks):
            fp_victim = random.choice(self.buckets[i])
            self.buckets[i].remove(fp_victim)
            self.buckets[i].append(fp)
            fp = fp_victim
            i = self._alt_index(i, fp)

            if len(self.buckets[i]) < self.bucket_size:
                self.buckets[i].append(fp)
                self.size += 1
                return True

        # Filter is considered full
        return False

    def __contains__(self, item):
        fp = self._fingerprint(item)
        i1 = self._hash(item)
        i2 = self._alt_index(i1, fp)
        return fp in self.buckets[i1] or fp in self.buckets[i2]

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

def main():
    # Total bits allocated for each filter
    TOTAL_BITS = 64 * 1024 * 1024 * 8  # 64 MB in bits
    print(f"Total bits allocated for each filter: {TOTAL_BITS}")

    # Target capacity: number of items
    capacity = 33_554_432

    # Bloom Filter Config
    bf_num_bits = TOTAL_BITS
    bf_bits_per_item = bf_num_bits / capacity
    bf_num_hashes = max(1, int(bf_bits_per_item * math.log(2)))
    print(f"\n--- Bloom Filter Configuration ---")
    print(f"Number of bits: {bf_num_bits}")
    print(f"Number of hashes: {bf_num_hashes}")
    print(f"Capacity: {capacity} items")
    print(f"Bits per item: {bf_bits_per_item:.2f} bits")

    bf = BloomFilter(num_bits=bf_num_bits, num_hashes=bf_num_hashes)

    # Cuckoo Filter Config
    cf_fingerprint_size = 16  # bits
    cf_bucket_size = 4
    cf_num_buckets = TOTAL_BITS // (cf_bucket_size * cf_fingerprint_size)
    cf_capacity = cf_num_buckets * cf_bucket_size
    if cf_capacity > capacity:
        cf_capacity = capacity
    else:
        capacity = cf_capacity  # Adjust capacity to match Cuckoo Filter
    cf_bits_per_item = TOTAL_BITS / cf_capacity
    print(f"\n--- Cuckoo Filter Configuration ---")
    print(f"Number of buckets: {cf_num_buckets}")
    print(f"Bucket size: {cf_bucket_size}")
    print(f"Fingerprint size: {cf_fingerprint_size} bits")
    print(f"Capacity: {cf_capacity} items")
    print(f"Bits per item: {cf_bits_per_item:.2f} bits")

    cf = CuckooFilter(num_buckets=cf_num_buckets, bucket_size=cf_bucket_size, fingerprint_size=cf_fingerprint_size)

    # Generate data for insertion
    num_elements = capacity
    print(f"\nGenerating {num_elements} random elements for insertion...")
    data = generate_random_strings(num_elements)

    # Insert elements into Bloom Filter
    print("\nInserting into Bloom Filter...")
    start_time = time.time()
    for item in data:
        bf.add(item)
    bf_insertion_time = time.time() - start_time
    print(f"Bloom Filter inserted {bf.size} elements in {bf_insertion_time:.2f} seconds.")

    # Insert elements into Cuckoo Filter
    print("\nInserting into Cuckoo Filter...")
    start_time = time.time()
    cf_inserted = 0
    for item in data:
        success = cf.add(item)
        if success:
            cf_inserted += 1
        else:
            break
    cf_insertion_time = time.time() - start_time
    print(f"Cuckoo Filter inserted {cf_inserted} elements in {cf_insertion_time:.2f} seconds.")

    # Prepare test data for false positive rate measurement
    num_test_items = 100000
    test_data = generate_random_strings(num_test_items)



    print("\nMeasuring False Positive Rates...")
    bf_fpr = measure_false_positive_rate(bf, "Bloom Filter",test_data , num_test_items)
    cf_fpr = measure_false_positive_rate(cf, "Cuckoo Filter",test_data, num_test_items)

    print("\n--- Summary ---")
    print(f"Bloom Filter:")
    print(f"  Inserted {bf.size} elements")
    print(f"  Bits per item: {bf_bits_per_item:.2f}")
    print(f"  False Positive Rate: {bf_fpr * 100:.4f}%")
    print(f"  Construction Time: {bf_insertion_time:.2f} seconds")
    print(f"Cuckoo Filter:")
    print(f"  Inserted {cf_inserted} elements")
    print(f"  Bits per item: {cf_bits_per_item:.2f}")
    print(f"  False Positive Rate: {cf_fpr * 100:.4f}%")
    print(f"  Construction Time: {cf_insertion_time:.2f} seconds")

if __name__ == "__main__":
    main()
