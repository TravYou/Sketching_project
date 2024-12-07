import time
import random
import string
import math
import mmh3

class QuotientFilter:
    def __init__(self, size):
        self.size = size
        self.table = [None] * size
        self.is_occupied = [False] * size
        self.is_continuation = [False] * size
        self.is_shifted = [False] * size

    def _hash(self, value):
        hash_value = hash(value)
        quotient = hash_value % self.size
        remainder = hash_value // self.size
        return quotient, remainder

    def insert(self, value):
        q, r = self._hash(value)
        if not self.is_occupied[q]:
            self.table[q] = r
            self.is_occupied[q] = True
        else:
            idx = q
            while self.is_occupied[idx]:
                idx = (idx + 1) % self.size
            self.table[idx] = r
            self.is_occupied[idx] = True
            self.is_shifted[idx] = True
            if idx != q:
                self.is_continuation[idx] = True

    def lookup(self, value):
        q, r = self._hash(value)
        idx = q
        while True:
            if not self.is_occupied[idx]:
                return False
            if self.table[idx] == r:
                return True
            if not self.is_continuation[idx]:
                break
            idx = (idx + 1) % self.size
        return False
    def __contains__(self, item):
        return self.lookup(item)

    def __str__(self):
        return str([
            (self.table[i], self.is_occupied[i], self.is_continuation[i], self.is_shifted[i])
            for i in range(self.size)
        ])
    
def generate_random_strings(num_strings, string_length=10):
    random_strings = set()
    while len(random_strings) < num_strings:
        s = ''.join(random.choices(string.ascii_letters + string.digits, k=string_length))
        random_strings.add(s)
    return list(random_strings)


def measure_false_positive_rate(filter_obj, filter_name, test_data, num_test_items):
    false_positives = 0
    for item in test_data:
        if item in filter_obj:
            false_positives += 1
    false_positive_rate = false_positives / num_test_items
    print(f"{filter_name} False Positive Rate: {false_positive_rate * 100:.4f}%")
    return false_positive_rate


def benchmark_quotient_filter():
    TOTAL_BITS = 8 * 1024 * 1024 * 8  #8MB in bits
    capacity = 1_000_000
    quotient_filter_size = capacity
    print(f"Total bits allocated for Quotient Filter: {TOTAL_BITS}")


    qf_num_slots = TOTAL_BITS // 32
    qf_capacity = qf_num_slots 
    qf_bits_per_item = TOTAL_BITS / qf_capacity

    print(f"\n--- Quotient Filter Configuration ---")
    print(f"Number of slots: {qf_num_slots}")
    print(f"Capacity: {qf_capacity} items")
    print(f"Bits per item: {qf_bits_per_item:.2f} bits")
    
    print(f"Generating {capacity} random elements for insertion")
    data = generate_random_strings(capacity)

    #initializing filter
    qf = QuotientFilter(size=quotient_filter_size)

    #inserting elements
    print("\nInserting into Quotient Filter...")
    start_time = time.time()
    for item in data:
        qf.insert(item)
    qf_insertion_time = time.time() - start_time
    print(f"Quotient Filter inserted {capacity} elements in {qf_insertion_time:.2f} seconds.")

    num_test_items = 10_000
    test_data = generate_random_strings(num_test_items)

    print("\nMeasuring False Positive Rates...")
    qf_fpr = measure_false_positive_rate(qf, "Quotient Filter", test_data, num_test_items)

    #print summary for table
    print("\n--- Summary for Quotient Filter ---")
    print(f"Inserted {capacity} elements")
    print(f"  Bits per item: {qf_bits_per_item:.2f}")
    print(f"False Positive Rate: {qf_fpr * 100:.4f}%")
    print(f"Construction Time: {qf_insertion_time:.2f} seconds")

if __name__ == "__main__":
    benchmark_quotient_filter()
