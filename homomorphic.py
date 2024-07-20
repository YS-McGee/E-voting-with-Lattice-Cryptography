import numpy as np

q = 3389  # Define the modulus

# Key Generation
def keygen(n):
    A_prime = np.random.randint(low=0, high=q, size=(n, n+1)) % q
    #print(A_prime)
    I_d = np.identity(n)
    #print(I_d)
    A = np.concatenate((A_prime, I_d), axis=1) % q
    #print(A)
    B = np.random.randint(low=0, high=q, size=(1, 2*n+1)) % q
    #print(B)
    C = np.concatenate((A, B), axis=0) % q
    #print(C)
    return C

# Commitment
def commit(C, m, n):
    r = np.random.randint(low=0, high=q, size=(2*n+1, 1)) % q
    m_vector = np.zeros((n+1, 1))
    m_vector[-1][0] = m
    c = (np.dot(C, r) + m_vector) % q
    return c, r

# Open Commitment
def open_commitment(C, c, r, m):
    m_prime_matrix = (c - np.dot(C, r)) % q
    m_prime = m_prime_matrix[-1, 0] % q
    if m_prime == m:
        return m_prime
    else:
        return False

# Sum Commitments
def sum_commitments(c1, c2):
    return (c1 + c2) % q

# Sum Random Vectors
def sum_random_vectors(r1, r2):
    return (r1 + r2) % q

# Example Usage
if __name__ == "__main__":
    n = 3  # dimension of matrix (d)
    C = keygen(n)

    # Commit to Number 1
    m1 = 3000  # Example number
    c1, r1 = commit(C, m1, n)

    # Commit to Number 2
    m2 = 400  # Example number
    c2, r2 = commit(C, m2, n)

    # Sum the Commitments
    c_sum = sum_commitments(c1, c2)
    r_sum = sum_random_vectors(r1, r2)
    m_sum = (m1 + m2) % q

    # Check commitments
    print(f"Commitment for number 1: {c1}")
    print(f"Commitment for number 2: {c2}")
    print(f"Sum of commitments: {c_sum}")

    # Open commitments
    open1 = open_commitment(C, c1, r1, m1)
    open2 = open_commitment(C, c2, r2, m2)
    open_sum = open_commitment(C, c_sum, r_sum, m_sum)

    print(f"Open commitment for number 1: {open1 }")
    print(f"Open commitment for number 2: {open2}")
    print(f"Open sum of commitments: {open_sum }")
