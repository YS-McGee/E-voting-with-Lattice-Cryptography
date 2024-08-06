import sys
import random
import numpy as np
# from sage.all import *

# For passing variables to cpp
import subprocess
import json

# Avoid printing scientific symbol, i.e. 2.869e+10
np.set_printoptions(suppress=True)

"""Commitment Constant Parameter"""
Nv = 5          # number of voters
N = 512         # Dimension of the lattice, change to 256 for real-life security
q = 1048583     # q > 2**20, such that f has an inverse in the ring
# q = 1073741827  # Large prime number for modulus as in Kyber (3389)


# Shamir's Secret Sharing implementation
def polynomial(coefs, x):
    """Evaluate a polynomial with given coefficients at x."""
    result = 0
    for coefficient in reversed(coefs):
        result = result * x + coefficient
    return result

def generate_shares(secret, num_shares, threshold):
    """Generate shares using Shamir's Secret Sharing scheme."""
    if threshold > num_shares:
        raise ValueError("Threshold cannot be greater than the number of shares")

    # coefs = [secret] + [random.randint(0, 100) for _ in range(threshold - 1)]
    coefs = [secret] + [np.random.randint(0, q+1) for _ in range(threshold - 1)]
    shares = [(i, polynomial(coefs, i)) for i in range(1, num_shares + 1)]
    print(f"secret: {secret}")
    print(f"coefs: {coefs}")
    print(f"shares; {shares}")
    return shares

# Lattice-based commitment scheme
def lattice_commitment(share, randomness, C):
    """Create a lattice-based commitment of the share and randomness."""
    # print(f"Share: {share}")
    # print(f"Randomness: {randomness}")
    # print(f"C: {C}")
    share_vector = np.zeros((N+1, 1))
    share_vector[-1][0] = share
    # print(f"share_vector: {share_vector}")
    commitment = (np.dot(C, randomness) + share_vector) % q
    # print(f"np.dot(C, randomness): {np.dot(C, randomness)}")
    return commitment

def open_commitment(trustees):
    """Open a lattice-based commitment to verify the share."""
    # opened = (np.dot(A, share) + randomness) % q
    # return np.array_equal(commitment, opened)
    for trustee in trustees:
        c_matrix = trustee['C_matrix']
        c_sum = trustee['sum_commitments']
        r_sum = trustee['sum_randomness']
        opened_share = (c_sum - np.dot(c_matrix, r_sum)) % q
        trustee['opened_share'] = opened_share[-1][0]

def voter(Nv, num_shares, threshold):
    voters = []

    for i in range(1, Nv + 1):
        # secret = random.randint(0, 1)
        secret = np.random.randint(0, 2)
        shares = generate_shares(secret, num_shares, threshold)
        
        voter_dict = {
            'voter_id': i,
            'secret_share': secret,
            'shares': shares
        }
        voters.append(voter_dict)

    return voters

def trustee(num_trustees):
    trustees = [{'trustee_id': i, 'received_shares': [], 'C_matrix': None, 'commitments': [], 'sum_commitments': None, 'sum_randomness': None, 'opened_share': None} for i in range(1, num_trustees + 1)]

    # Key Generation
    def keygen(N):
        # A_prime = random.randint(low=0, high=q, size=(N, N+1)) % q
        A_prime = np.random.randint(low=0, high=q+1, size=(N, N+1)) % q
        #print(A_prime)
        I_d = np.identity(N)
        #print(I_d)
        A = np.concatenate((A_prime, I_d), axis=1) % q
        #print(A)
        B = np.random.randint(low=0, high=q, size=(1, 2*N+1)) % q
        #print(B)
        C = np.concatenate((A, B), axis=0) % q                              # C's dimension is (n+1)X(2n+1)
        #print(C)
        return C

    # Insert secret matrix for each trustee
    for trustee in trustees:
        trustee['C_matrix'] = keygen(N)

    return trustees

def commit_shares(trustees, A):
    # print(f"\n[TRUSTEE] -- {trustees} \n")
    for trustee in trustees:
        # print(f"\n[TRUSTEE] -- {trustee}")
        # sum_commitment = np.zeros(A.shape[0])
        # # print(sum_commitment)
        # sum_randomness = np.zeros(A.shape[0])

        C = trustee['C_matrix']
        test_sum_commitment = np.zeros(1)
        test_sum_randomness = np.zeros(1)
        # print(test_randomness)

        # for share in trustee['received_shares']:
        #     share_value = share['share_value']
        #     randomness = np.random.randint(0, q, size=A.shape[0])  # Randomness for commitment
        #     commitment = lattice_commitment(share_value, randomness, A, q)
        #     trustee['commitments'].append({
        #         'voter_id': share['voter_id'],
        #         'commitment': commitment,
        #         'randomness': randomness
        #     })
        #     sum_commitment = (sum_commitment + commitment) % q
        #     sum_randomness = (sum_randomness + randomness) % q

        for share in trustee['received_shares']:
            share_value = share['share_value']
            # print(f"share_value: {share_value}")
            randomness = np.random.randint(low=0, high=q, size=(2*N+1, 2*N+1))  # Randomness for commitment
            #print(f"randomness:\n {randomness}")
            commitment = lattice_commitment(share_value, randomness, C)
            # print(f"commitment: {commitment}")
            trustee['commitments'].append({
                'voter_id': share['voter_id'],
                'commitment': commitment,
                'randomness': randomness
            })
            test_sum_commitment = (test_sum_commitment + commitment) % q
            test_sum_randomness = (test_sum_randomness + randomness) % q

        trustee['sum_commitments'] = test_sum_commitment
        trustee['sum_randomness'] = test_sum_randomness
        
        # print(f"\n[DEBUG] -- {trustee}")

def sum_commitments(trustees):
    for trustee in trustees:
        sum_commitment = np.zeros_like(trustee['commitments'][0]['commitment'])
        for commitment in trustee['commitments']:
            sum_commitment = (sum_commitment + commitment['commitment']) % q
        trustee['sum_commitments'] = sum_commitment
        # print(f"\n[TRUSTEE] --\n{trustee}")

# def open_sum(trustees, A, q):
#     for trustee in trustees:
#         #print(f"trustee {trustee}")
#         sum_commitment = trustee['sum_commitments']
#         sum_randomness = trustee['sum_randomness']
#         sum_share = sum([share['share_value'] for share in trustee['received_shares']]) % q
#         if open_commitment(sum_commitment, sum_share, sum_randomness, A, q):
#             # print(f"Trustee {trustee['trustee_id']} successfully opened the sum commitment.")
#             pass
#         else:
#             print(f"Trustee {trustee['trustee_id']} failed to open the sum commitment.")
#     total_sum_share = sum([sum([share['share_value'] for share in trustee['received_shares']]) for trustee in trustees]) % q
#     total_commitments_sum = sum([trustee['sum_commitments'] for trustee in trustees]) % q
#     print(f"Total sum of commitments from trustees:\n {total_commitments_sum}")
#     print(f"Total sum of actual votes: {total_sum_share}")

def distribute_shares(voters, trustees):
    for voter in voters:
        for share in voter['shares']:
            share_id, share_value = share
            trustees[share_id - 1]['received_shares'].append({
                'voter_id': voter['voter_id'],
                'share_value': share_value
            })

# Lagrange interpolation function to compute sum of votes
def lagrange_interpolation(shares, x):
    print(f"shares: {shares}")
    total = 0
    n = len(shares)
    for i in range(n):
        xi, yi = shares[i]
        term = yi
        for j in range(n):
            if i != j:
                xj, _ = shares[j]
                # Compute the modular inverse
                term *= (x - xj) * pow(xi - xj, -1, q)
                term %= q
        total += term
        total %= q
    return total

# Using SageMath
# def lagrange_interpolation(shares, x, q):
#     total = 0
#     n = len(shares)
#     x = Integer(x)
#     q = Integer(q)
    
#     for i in range(n):
#         xi, yi = Integer(shares[i][0]), Integer(shares[i][1])
#         term = yi
#         for j in range(n):
#             if i != j:
#                 xj, _ = Integer(shares[j][0]), Integer(shares[j][1])
#                 # Compute the modular inverse using SageMath's .inverse_mod() method
#                 try:
#                     inv = (xi - xj).inverse_mod(q)
#                     term *= (x - xj) * inv
#                     term %= q
#                 except ZeroDivisionError:
#                     print(f"No modular inverse for {xi} - {xj} modulo {q}.")
#                     return None
#         total += term
#         total %= q
#     return total

def main():
    try:
        if len(sys.argv) > 2:
            print(f"len(sys.argv) = {len(sys.argv)}")
            print("Invalid input. Please enter a valid integer.")
            return False
        elif len(sys.argv) == 2:
            user_input = sys.argv[1]
            try:
                number = int(user_input)
                if number < 1 or number > 10:
                    print("Invalid input. Please enter a number between 1 and 10.")
                    return
                num_voters = int(input("Enter the number of voters: "))

                ######################## Nv = int(input("Enter the number of voters: ")) ########################

                print(f"Number of voters: {num_voters}")
            except ValueError:
                print("Invalid input. Please enter a valid integer.")
                return
        elif len(sys.argv) > 1:
            pass                                # No input provided, just pass
    except ValueError:
        print("Invalid input. Please enter a valid integer.")


    threshold = 2                               # minimum number of available trustees
    num_shares = threshold * 2 - 1

    num_trustees = num_shares                   # Number of trustees should be equal to num_shares
    print(f"Number of Trustees: {num_trustees}")

    # Lattice parameters
    A = np.random.randint(0, 100, size=(N, 1))  # Generating a random matrix A with elements in [0, 100)

    voters_list = voter(Nv, num_shares, threshold)
    trustees_list = trustee(num_trustees)
    # print(f"trustees_list {trustees_list}")
    
    distribute_shares(voters_list, trustees_list)
    # print(f"trustees_list {trustees_list}")

    # Trustees Main Operations
    commit_shares(trustees_list, A)
    #print(f"[DEBUG] -- {trustees_list}")

    # Each trsutee open the sum of shares
    open_commitment(trustees_list)

    # sum_commitments(trustees_list, q)
    # open_sum(trustees_list, A, q)

    for t in trustees_list:
        print(f"[trustees_list] -- {t['trustee_id']}, {t['received_shares']}")

    for v in voters_list:
        print(f"[voters_list] -- {v}")
    # print(f"voters_list: {voters_list}")

    trustee_shares = [(t['trustee_id'], t['opened_share']) for t in trustees_list]
    shares = [(int(x), int(y)) for x, y in trustee_shares]
    secret_sum = lagrange_interpolation(shares, 0)
    print("Reconstructed Secret Sum:", secret_sum)
    # print(f"type trustees_list {type(trustees_list)}")
    # print(f"type voters_list {type(voters_list)}")

    # arg1 = "Hello"
    # arg2 = 42
    # subprocess.run(["./Lattice-IBE-master/IBE", arg1, str(arg2)])

    # data = voters_list

    data_str = json.dumps(voters_list)
    subprocess.run(["./Lattice-IBE-master/IBE", data_str])

if __name__ == "__main__":
    main()
