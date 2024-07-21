import sys
import random
import numpy as np

np.set_printoptions(suppress=True)

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

    coefs = [secret] + [random.randint(0, 100) for _ in range(threshold - 1)]
    shares = [(i, polynomial(coefs, i)) for i in range(1, num_shares + 1)]
    return shares

# Lattice-based commitment scheme
def lattice_commitment(share, randomness, C, n, q):
    """Create a lattice-based commitment of the share and randomness."""
    #print(f"Share: {share}")
    # print(f"Randomness: {randomness}")
    # print(f"C: {C}")
    share_vector = np.zeros((n+1, 1))
    share_vector[-1][0] = share
    commitment = (np.dot(C, randomness) + share_vector) % q

    return commitment

def open_commitment(commitment, share, randomness, A, q):
    """Open a lattice-based commitment to verify the share."""
    opened = (np.dot(A, share) + randomness) % q
    return np.array_equal(commitment, opened)

def voter(Nv, num_shares, threshold):
    voters = []

    for i in range(1, Nv + 1):
        secret = random.randint(0, 1)
        shares = generate_shares(secret, num_shares, threshold)
        
        voter_dict = {
            'voter_id': i,
            'secret_share': secret,
            'shares': shares
        }
        voters.append(voter_dict)

    return voters

def trustee(num_trustees, n, q):
    trustees = [{'trustee_id': i, 'received_shares': [], 'C_matrix': None, 'commitments': [], 'sum_commitments': None, 'sum_randomness': None} for i in range(1, num_trustees + 1)]

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
        C = np.concatenate((A, B), axis=0) % q                              # C's dimension is (d+1)X(2d+1)
        #print(C)
        return C

    # Insert secret matrix for each trustee
    for trustee in trustees:
        trustee['C_matrix'] = keygen(n)

    return trustees

def commit_shares(trustees, A, n, q):
    # print(f"\n[TRUSTEE] -- {trustees} \n")
    for trustee in trustees:
        # print(f"\n[TRUSTEE] -- {trustee}")
        sum_commitment = np.zeros(A.shape[0])
        # print(sum_commitment)
        sum_randomness = np.zeros(A.shape[0])

        C = trustee['C_matrix']
        test_sum_commitment = np.zeros(C.shape[1])
        test_sum_randomness = np.zeros(C.shape[1])
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
            randomness = np.random.randint(low=0, high=q, size=(2*n+1, 1))  # Randomness for commitment
            print(f"randomness:\n {randomness}")
            commitment = lattice_commitment(share_value, randomness, C, n, q)
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
        
        print(f"\n[TRUSTEE] -- {trustee}")

def sum_commitments(trustees, q):
    for trustee in trustees:
        sum_commitment = np.zeros_like(trustee['commitments'][0]['commitment'])
        for commitment in trustee['commitments']:
            sum_commitment = (sum_commitment + commitment['commitment']) % q
        trustee['sum_commitments'] = sum_commitment

def open_sum(trustees, A, q):
    for trustee in trustees:
        print(f"trustee {trustee}")
        sum_commitment = trustee['sum_commitments']
        sum_randomness = trustee['sum_randomness']
        sum_share = sum([share['share_value'] for share in trustee['received_shares']]) % q
        if open_commitment(sum_commitment, sum_share, sum_randomness, A, q):
            # print(f"Trustee {trustee['trustee_id']} successfully opened the sum commitment.")
            pass
        else:
            print(f"Trustee {trustee['trustee_id']} failed to open the sum commitment.")
    total_sum_share = sum([sum([share['share_value'] for share in trustee['received_shares']]) for trustee in trustees]) % q
    total_commitments_sum = sum([trustee['sum_commitments'] for trustee in trustees]) % q
    print(f"Total sum of commitments from trustees:\n {total_commitments_sum}")
    print(f"Total sum of actual votes: {total_sum_share}")

def distribute_shares(voters, trustees):
    for voter in voters:
        for share in voter['shares']:
            share_id, share_value = share
            trustees[share_id - 1]['received_shares'].append({
                'voter_id': voter['voter_id'],
                'share_value': share_value
            })

def main():

    Nv = 3          # number of voters
    
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

    # Lattice parameters
    n = 3                                       # Dimension of the lattice, change to 256 for real-life security
    A = np.random.randint(0, 100, size=(n, 1))  # Generating a random matrix A with elements in [0, 100)
    q = 3389                                    # Large prime number for modulus as in Kyber

    voters_list = voter(Nv, num_shares, threshold)
    for v in voters_list:
        print(f"voters_list {v}")
    
    trustees_list = trustee(num_trustees, n, q)
    # print(f"trustees_list {trustees_list}")
    
    distribute_shares(voters_list, trustees_list)
    # print(f"trustees_list {trustees_list}")

    # Trustees Main Operations
    commit_shares(trustees_list, A, n, q)
    sum_commitments(trustees_list, q)
    open_sum(trustees_list, A, q)

    for v in voters_list:
        print(v)

    for t in trustees_list:
        print(t)

if __name__ == "__main__":
    main()
