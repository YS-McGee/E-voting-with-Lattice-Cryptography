import numpy as np
from numpy.fft import fft, ifft
from scipy.stats import norm as scipy_norm
from sympy import symbols, ZZ, Poly
import sys
from sage.all import *

# Constants
N = 512
q = 2**30
sigma_1 = 0.84932180028801904272150283410288961971514109378435394286159953238339383120795466719298223538163406787061691601172910413284884326532697308797136114023
LDRMX = 2**31 - 1
log_2 = np.log(2)
omega = np.exp(2j * np.pi / N)
omega_1 = np.exp(-2j * np.pi / N)

# For Poly
# x = symbols('x')

# for Sagemath
R = PolynomialRing(ZZ, 'x')
Rx = R.gen()

def sample0(alpha):
    """Samples from distribution D_{sigma_2}^+."""
    if (alpha & 1) == 0:
        return 0
    i = 1
    k = 1
    mask = 0
    while i < 1000:
        aux = alpha & mask
        alpha = alpha >> k
        if aux:
            return sample0(alpha)
        else:
            if (alpha & 1) == 0:
                return i
        i += 1
        k += 2
        mask = (mask << 2) | 6
    print("ERROR")
    return 999999

def sample1(k):
    """Samples from distribution D_{k*sigma_2}^+."""
    alpha = np.random.randint(0, LDRMX)
    x = sample0(alpha)
    y = np.random.randint(0, k)
    z = k * x + y
    w = y * ((z << 1) - y)
    bravo = LDRMX / np.exp(w * log_2 / (k * k))
    alpha = np.random.randint(0, LDRMX)
    if alpha > bravo:
        return sample1(k)
    else:
        return z

def sample2(k):
    """Samples from distribution D_{k*sigma_2}."""
    while True:
        alpha = np.random.randint(0, LDRMX)
        x = sample1(k)
        if x != 0 or (alpha & 1) == 1:
            alpha >>= 1
            signe = 1 - 2 * (alpha & 1)
            return x * signe

def sample3(sigma):
    """Samples from distribution D_{sigma}."""
    k = int(np.ceil(sigma / sigma_1))
    while True:
        x = sample2(k)
        alpha = np.random.rand()  # Equivalent to ((double)rand()) / LDRMX in C++
        bravo = np.exp(-x * x * (1 / (2 * sigma * sigma) - 1 / (2 * k * k * sigma_1 * sigma_1)))
        assert bravo <= 1
        if alpha < bravo:
            return x

def sample4(c, sigma):
    """Samples from distribution D_{c,sigma}."""
    intc = int(np.floor(c))
    fracc = c - intc
    denom = 1 / (2 * sigma * sigma)
    while True:
        x = sample3(sigma)
        flip = np.random.randint(0, 2)
        x += flip
        bravo = np.exp(-(x - fracc) * (x - fracc) * denom) / (np.exp(-x * x * denom) + np.exp(-(x - 1) * (x - 1) * denom))
        assert bravo < 1
        alpha = np.random.rand()
        if alpha < bravo:
            return x + intc

def sample_polynomial(N, sigma):
    """Sample a polynomial with coefficients from a discrete Gaussian distribution."""
    # return [sample_discrete_gaussian(sigma) for _ in range(N)]
    coeffs = np.array([sample3(sigma) for _ in range(N)], dtype=np.int64)
    coeffs[-1] |= 1  # Ensure the leading coefficient is non-zero
    return coeffs

def FFTStep(f, N, w0):
    """Recursive FFT step."""
    if N == 1:
        return f
    elif N == 2:
        return np.array([f[0] + 1j * f[1], f[0] - 1j * f[1]])
    else:
        assert N % 2 == 0
        f0 = f[0::2]
        f1 = f[1::2]

        w02 = w0 * w0
        wk = w0
        f0_fft = FFTStep(f0, N // 2, w02)
        f1_fft = FFTStep(f1, N // 2, w02)
        f_fft = np.zeros(N, dtype=complex)
        for k in range(N):
            f_fft[k] = f0_fft[k % (N // 2)] + wk * f1_fft[k % (N // 2)]
            wk *= w02
        return f_fft

def ZZXToFFT(f):
    """Convert polynomial to its FFT representation."""
    f_double = np.array(f, dtype=float)
    assert len(f) == N
    return FFTStep(f_double, N, omega)

def FFTRealReverse(f_fft):
    """Perform inverse FFT and return real part of the result."""
    fprime = np.fft.ifft(f_fft)
    return fprime.real

def gram_schmidt_norm(f, g):
    # Norm of (g, -f)
    norm_1 = np.sqrt(np.sum(f**2 + g**2))
    # print(f"norm_1: {norm_1}")

    f_fft = ZZXToFFT(f)
    g_fft = ZZXToFFT(g)

    F = np.zeros(N, dtype=complex)
    G = np.zeros(N, dtype=complex)

    for i in range(N):
        denominator = f_fft[i]*f_fft[N-1-i] + g_fft[i]*g_fft[N-1-i]
        F[i] = f[i] / denominator
        G[i] = g[i] / denominator

    Foxtrot = FFTRealReverse(F)
    Golf = FFTRealReverse(G)

    norm_2 = q * np.sqrt(np.sum(Foxtrot**2 + Golf**2))
    # print(f"norm_2: {norm_2}")

    return max(norm_1, norm_2)

def Cyclo():
    """Generate the cyclotomic polynomial."""
    return Rx**N + 1

def numpy_to_sage_poly(np_array):
    """Convert a numpy array to a SageMath polynomial."""
    coeffs = np_array.tolist()
    return R(coeffs[:])

def sage_to_numpy(sage_poly):
    # Extract coefficients and reverse to match numpy polynomial coefficient order (highest degree first)
    coefficients = sage_poly.coefficients(sparse=False)
    
    # Convert to numpy array
    numpy_array = np.array(coefficients)
    
    return numpy_array

def reduce_poly_mod(poly):
    return [(poly[i] - sum(poly[j] for j in range(N, len(poly)) if j - i == N)) for i in range(N)]

def poly_mul_mod(a, b, mod_poly):
    """Multiply two polynomials and reduce modulo another polynomial."""
    return (a * b) % mod_poly

def pair_gcd(f, g):
    """SAGE env. Compute the GCD of f and g with respect to the cyclotomic polynomial phi using extended GCD."""
    # print(f"f type: {type(f)}")
    sage_f = numpy_to_sage_poly(f)
    sage_g = numpy_to_sage_poly(g)

    phi = Cyclo()

    """
    First compute Rf and test GCD(Rf,q)
    """
    # Compute extended GCD
    Rf, rho_f, _ = xgcd(sage_f, phi)
    if gcd(Rf, q) != 1:
        return False, None, None, None, None

    """
    Compute Rg and test GCD(Rf,Rg)
    """
    Rg, rho_g, _ = xgcd(sage_g, phi)
    gcd_Rf_Rg, u, v = xgcd(Rf, Rg)

    # Check GCD(Rf,Rg)
    if gcd_Rf_Rg != 1:
        return False, None, None, None, None
    
    return True, u, v, rho_f, rho_g

def unique_reverse(coeffs):
    """Compute the unique polynomial f_bar."""
    return np.array([coeffs[0]] + [-coeffs[N-i] for i in range(1, N)])

def poly_mult_mod(f, g):
    # polynomial ring multiplication 
    poly = [sum(f[i] * g[j] for i in range(len(f)) for j in range(len(g)) if i + j == k) for k in range(len(f) + len(g) - 1)]
    return reduce_poly_mod(poly)

def poly_mult_fft(f, g):
    size = len(f) + len(g) - 1
    size_padded = 2 ** int(np.ceil(np.log2(size)))  # Next power of 2
    f_padded = np.pad(f, (0, size_padded - len(f)))
    g_padded = np.pad(g, (0, size_padded - len(g)))
    
    # FFT
    fft_f = np.fft.fft(f_padded)
    fft_g = np.fft.fft(g_padded)
    
    # Pointwise multiplication
    fft_result = fft_f * fft_g
    
    # Inverse FFT
    result = np.fft.ifft(fft_result)
    result = np.round(result).astype(int)

    # Polynomial Ring
    # return reduce_poly_mod(result[:size])
    return result[:size]

def multiply_large_polynomials(poly1, poly2):
    """
    Multiply two polynomials with extremely large coefficients using SageMath's arbitrary-precision integers.
    
    Args:
        poly1 (list): Coefficients of the first polynomial.
        poly2 (list): Coefficients of the second polynomial.
    
    Returns:
        list: Coefficients of the resulting polynomial after multiplication.
    """
    # Define a polynomial ring with ZZ coefficients
    R = PolynomialRing(ZZ, 'x')
    x = R.gen()
    
    # Convert lists to SageMath polynomials
    P1 = R(poly1)
    P2 = R(poly2)
    
    # Perform polynomial multiplication
    P_result = P1 * P2
    
    # Convert the result back to a list of coefficients
    result_coeffs = P_result.coefficients(sparse=False)
    
    return reduce_poly_mod(result_coeffs)

def compute_k(f, g, F, G):
    """Compute the reduction coefficient k."""
    f_bar = unique_reverse(f).tolist()
    g_bar = unique_reverse(g).tolist()

    f = f.tolist()
    g = g.tolist()

    alpha = multiply_large_polynomials(F, f_bar)
    bravo = multiply_large_polynomials(G, g_bar)
    charlie = multiply_large_polynomials(f, f_bar)
    delta = multiply_large_polynomials(g, g_bar)

    nume = [a + b for a, b in zip(alpha, bravo)]
    deno = [c + d for c, d in zip(charlie, delta)]
    # print(f"nume: {nume}")
    # print(f"deno: {deno}")

    # Convert to Sage polynomials for division and GCD operations
    # sage_nume = numpy_to_sage_poly(np.array(nume))
    sage_deno = numpy_to_sage_poly(np.array(deno))

    # Perform extended GCD to find gcd, inverse of deno mod phi
    phi = Cyclo()
    gcd, inv_deno, _ = xgcd(sage_deno, phi)

    inv_deno = sage_to_numpy(inv_deno).tolist()

    # # print(f"sage_to_numpy(inv_deno): {sage_to_numpy(inv_deno)}")
    k = multiply_large_polynomials(nume, inv_deno)
    k = numpy_to_sage_poly(np.array(k))
    # SageMath Opeartion
    k = k // gcd
    # print(f"k // gcd: {k}")

    # Ensure k has length N
    k_coeffs = k.list()
    if len(k_coeffs) < N:
        k_coeffs += [0] * (N - len(k_coeffs))
    k_coeffs = k_coeffs[:N]

    # print(f"k_coeffs: {k_coeffs}")

    return k_coeffs

def Inverse(f):
    
    # Convert f and phi to polynomials in R
    f = R(f)
    phi = Cyclo()
    
    # Compute the extended GCD of f and phi
    gcd, rho_f, iphi = xgcd(f, phi)
    inv_f = inverse_mod(gcd, q)

    # print(f"gcd: {gcd}")
    # print(f"inv_f: {inv_f}")
    
    # Return the inverse polynomial
    return inv_f * rho_f

def mpk_gen(f, g):
    # Define the polynomial ring over the finite field GF(q)
    R = PolynomialRing(GF(q), 'x')
    x = R.gen()

    print(f"f: {f}")
    print(f"g: {g}")

    f = R(f)
    g = R(g)

    f_inv = Inverse(f)
    f_inv = sage_to_numpy(f_inv).tolist()
    # print(f"f: {f}")
    # print(f"g: {g}")
    # print(f"R(g): {R(g)}")
    print(f"f_inv: {f_inv}")
    h = [a%q for a in multiply_large_polynomials(g, f_inv)]

    print(f"h: {h}")
    return h

def key_generation(N, q):
    sigma_f = 1.17 * np.sqrt(q / (2 * N))
    
    i = 1
    while True:
        # print(f"i: {i}")
        f = sample_polynomial(N, sigma_f)
        g = sample_polynomial(N, sigma_f)
        # print(f"f: {f}")
        # print(f"g: {g}")
        
        norm = gram_schmidt_norm(f, g)
        valid, u, v, sage_rho_f, sage_rho_g = pair_gcd(f, g)

        if norm <= 1.17 * np.sqrt(q) and valid:
            # print(f"norm: {norm}")
            print("norm and gcd ok!")
            break
        i += 1

    # Convert sage ring polynomial to numpy.ndarray
    rho_f = sage_to_numpy(sage_rho_f)
    rho_g = sage_to_numpy(sage_rho_g)

    # print(f"q: {q}")
    # print(f"v: {v}")
    # print(f"rho_g: {rho_g}")

    F = [-q*v*num for num in rho_g]
    G = [q*u*num for num in rho_f]

    k = compute_k(f, g, F, G)
    while R(k).degree() >= 0:
        f = f.tolist()
        g = g.tolist()

        F = [a-b for a, b in zip(F, multiply_large_polynomials(k, f))]
        G = [a-b for a, b in zip(G, multiply_large_polynomials(k, g))]

        f = np.array(f)
        g = np.array(g)

        k = compute_k(f, g, F, G)
        # print(f"k: {k}")
    
    f = f.tolist()
    g = g.tolist()

    q_test = [a-b for a, b in zip(multiply_large_polynomials(f, G), multiply_large_polynomials(g, F))]
    q_test = q_test[0]
    if q == q_test:
        print(f"KeyGen Successful! f*G - g*F = q = {q_test}")
    else:
        print("FAILED f*G - g*F != q")
        sys.exit()

    MPK = mpk_gen(f, g)

    print(f"MPK: {MPK}")
    
    return 1

if __name__ == "__main__":
    # Generate key
    valid = key_generation(N, q)


    

    