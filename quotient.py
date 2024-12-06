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
            # Insert the new value
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

    def __str__(self):
        return str([
            (self.table[i], self.is_occupied[i], self.is_continuation[i], self.is_shifted[i])
            for i in range(self.size)
        ])


if __name__ == "__main__":
    qf = QuotientFilter(size=16)

    qf.insert("apple")
    qf.insert("banana")
    qf.insert("cherry")

    
    print("apple:", qf.lookup("apple"))  # True
    print("grape:", qf.lookup("grape"))  # False



