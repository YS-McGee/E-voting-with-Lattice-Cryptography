import numpy as np
from scipy.stats import norm

def create_cdt(sigma, length):
    """Create a CDT for the discrete Gaussian distribution."""
    x = np.arange(0, length)
    cdf = norm.cdf(x + 0.5, scale=sigma) - norm.cdf(x - 0.5, scale=sigma)
    cdf[1:] = 2 * cdf[1:]  # Double the probability for positive x
    cdt = np.cumsum(cdf)
    return cdt / cdt[-1]  # Normalize to create the CDT

def sample_cdt(cdt):
    """Sample from the CDT using binary search."""
    u = np.random.random()
    index = np.searchsorted(cdt, u)
    return index

def sample_discrete_gaussian(sigma, length=100):
    """Sample from a discrete Gaussian distribution using CDT."""
    cdt = create_cdt(sigma, length)
    x = sample_cdt(cdt)
    return x if np.random.random() < 0.5 else -x

def sample_polynomial(N, sigma):
    """Sample a polynomial with coefficients from a discrete Gaussian distribution."""
    return [sample_discrete_gaussian(sigma) for _ in range(N)]

def gram_schmidt_norm(B):
    """Compute the Gram-Schmidt norm of B."""
    B = B.astype(np.float64)  # Ensure B is of type float64
    B_star = np.zeros_like(B)
    for i in range(B.shape[1]):
        B_star[:, i] = B[:, i]
        for j in range(i):
            B_star[:, i] -= np.dot(B[:, i], B_star[:, j]) / np.dot(B_star[:, j], B_star[:, j]) * B_star[:, j]
    return max(np.linalg.norm(B_star[:, i]) for i in range(B.shape[1]))

def compute_norm(f, g, q):
    """Compute the norm condition as described in the algorithm."""
    f = np.array(f)
    g = np.array(g)
    f_star = f + g * f
    g_star = g + f * g
    B = np.vstack((f, g))
    B_star = np.vstack((q * f, f_star, q * g, g_star))
    return max(np.linalg.norm(B[i]) for i in range(B.shape[0]))

def key_generation(N, q):
    sigma_f = 1.17 * np.sqrt(q / (2 * N))
    
    while True:
        f = sample_polynomial(N, sigma_f)
        g = sample_polynomial(N, sigma_f)
        
        norm_val = compute_norm(f, g, q)
        if norm_val <= 1.17 * np.sqrt(q):
            break

    return f, g

if __name__ == "__main__":
    # Example parameters
    N = 4
    q = 2**23

    # Generate key
    f, g = key_generation(N, q)
    print("f:", f)
    print("g:", g)
