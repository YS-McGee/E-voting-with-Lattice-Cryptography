import numpy as np

# Votes
votes = [1, 0, 1, 1, 0]

# Modulus
q = 3389

# Generate shares for each vote
def generate_shares(secret, num_shares, threshold, q):
    coefficients = np.random.randint(0, q, threshold - 1).tolist()
    coefficients.insert(0, secret)
    shares = [(i, sum([coeff * (i ** idx) for idx, coeff in enumerate(coefficients)]) % q) for i in range(1, num_shares + 1)]
    return shares

# Generate shares for all votes
all_shares = [generate_shares(vote, 3, 2, q) for vote in votes]
print("All Shares:", all_shares)

# Sum shares for each trustee
summed_shares = [sum([all_shares[j][i][1] for j in range(len(all_shares))]) % q for i in range(len(all_shares[0]))]
print("Summed Shares:", summed_shares)

# Define the Lagrange interpolation function
def lagrange_interpolation(shares, x, q):
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

# Reconstruct the secret
trustee_shares = [(i + 1, summed_shares[i]) for i in range(len(summed_shares))]
secret_sum = lagrange_interpolation(trustee_shares, 0, q)
print("Reconstructed Secret Sum:", secret_sum)
