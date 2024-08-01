import numpy as np
from numpy.fft import fft, ifft
from scipy.stats import norm as scipy_norm
from sympy import symbols, ZZ, Poly
# from egcd import egcd
from sage.all import *
from sage.arith.misc import *

# Constants
N = 4
q = 2**30
sigma_1 = 0.84932180028801904272150283410288961971514109378435394286159953238339383120795466719298223538163406787061691601172910413284884326532697308797136114023
LDRMX = 2**31 - 1
log_2 = np.log(2)
omega = np.exp(2j * np.pi / N)
omega_1 = np.exp(-2j * np.pi / N)

# For Poly
x = symbols('x')

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
    
# def Cyclo():
#     """Generate the cyclotomic polynomial."""
#     coeffs = [1] + [0]*(N-1) + [1]
#     phi = Poly(coeffs, x, domain=ZZ)
#     return phi

# def pair_gcd(f, g):
#     """Compute the GCD of f and g with respect to the cyclotomic polynomial phi using egcd."""
#     phi = Cyclo()

#     # Convert numpy arrays to Poly objects
#     f = Poly(f[::-1], x, domain=ZZ)
#     g = Poly(g[::-1], x, domain=ZZ)
#     print(f"f: {f}")
#     print(f"g: {g}")

#     # Transform polynomials to their coefficient arrays
#     f_coeffs = np.array(f.all_coeffs()[::-1], dtype=np.int64)
#     g_coeffs = np.array(g.all_coeffs()[::-1], dtype=np.int64)
#     print(f"f_coeffs: {f_coeffs}")
#     print(f"g_coeffs: {g_coeffs}")
#     phi_coeffs = np.array(phi.all_coeffs()[::-1], dtype=np.int64)

#     # Compute the GCD using egcd for each coefficient
#     res_f_gcd = np.zeros_like(f_coeffs)
#     for i in range(len(f_coeffs)):
#         res_f_gcd[i], _, _ = egcd(int(f_coeffs[i]), int(phi_coeffs[i]))
#     print(f"res_f_gcd: {res_f_gcd}")
#     if np.gcd.reduce(res_f_gcd) != 1:
#         pgcd = ZZ(0)
#         alpha = None
#         beta = None
#         rho_f = None
#         rho_g = None
#     else:
#         res_g_gcd = np.zeros_like(g_coeffs)
#         for i in range(len(g_coeffs)):
#             res_g_gcd[i], _, _ = egcd(int(g_coeffs[i]), int(phi_coeffs[i]))
#         pgcd = np.gcd.reduce(res_f_gcd)
#         alpha = ZZ(1)
#         beta = ZZ(1)
#         rho_f = Poly(res_f_gcd[::-1], x, domain=ZZ)
#         rho_g = Poly(res_g_gcd[::-1], x, domain=ZZ)

#     return pgcd, alpha, beta, rho_f, rho_g

def Cyclo():
    """Generate the cyclotomic polynomial."""
    return Rx**N + 1

def numpy_to_sage_poly(np_array):
    """Convert a numpy array to a SageMath polynomial."""
    coeffs = np_array.tolist()
    return R(coeffs[::-1])

def reduce_poly_mod(poly):
    return [(poly[i] - sum(poly[j] for j in range(N, len(poly)) if j - i == N)) for i in range(N)]

def poly_mul_mod(a, b, mod_poly):
    """Multiply two polynomials and reduce modulo another polynomial."""
    return (a * b) % mod_poly

def pair_gcd(f, g):
    """Compute the GCD of f and g with respect to the cyclotomic polynomial phi using extended GCD."""
    f = numpy_to_sage_poly(f)
    g = numpy_to_sage_poly(g)

    phi = Cyclo()

    """
    First compute Rf and test GCD(Rf,q)
    """
    # Compute extended GCD
    gcd_f, rho_f, _ = xgcd(f, phi)
    Rf = poly_mul_mod(rho_f, f, phi)
    # Reduce to modulo q
    Rf = Rf.constant_coefficient() % q

    if gcd(Rf, q) != 1:
        return False, None, None, None, None

    """
    Compute Rg and test GCD(Rf,Rg)
    """
    gcd_g, rho_g, _ = xgcd(g, phi)
    Rg = poly_mul_mod(rho_g, g, phi)
    # Reduce to modulo q
    Rg = Rg.constant_coefficient() % q

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

def compute_k(f, g, F, G):
    """Compute the reduction coefficient k."""
    # print(f"f: {f}")
    # print(f"g: {g}")
    # v =poly_mult_mod(f, g)
    # print(f"f*g: {v}")
    # print(f"mod f*g: {reduce_poly_mod(v)}")
    
    f_bar = unique_reverse(f)
    g_bar = unique_reverse(g)
    print(f"f_bar type: {type(f_bar)}")
    print(f"F type: {type(F)}")

    nume = reduce_poly_mod(poly_mult_mod(F, f_bar) + poly_mult_mod(G, g_bar))
    print(f"nume: {nume}")

    return 1

    # R = PolynomialRing(ZZ, 'x')
    # x = R.gen()

    # phi = PolynomialRing(ZZ, 'x').gen()**N + 1  # cyclotomic polynomial

    # # Convert numpy arrays to Sage polynomials
    # f_poly = R(list(f[::-1]))
    # g_poly = R(list(g[::-1]))

    #  # Compute f_bar and g_bar
    # f_coeffs = np.array(f)
    # g_coeffs = np.array(g)
    # f_bar_coeffs = unique_reverse(f_coeffs)
    # g_bar_coeffs = unique_reverse(g_coeffs)
    # f_bar_poly = R(list(f_bar_coeffs[::-1]))
    # g_bar_poly = R(list(g_bar_coeffs[::-1]))

    # # Compute numerator and denominator
    # numerator = reduce_poly_mod(f_bar_poly * F + g_bar_poly * G)
    # denominator = reduce_poly_mod(f_poly * f_bar_poly + g_poly * g_bar_poly)

    # # Compute the inverse of the denominator modulo phi
    # a, iden, _ = xgcd(denominator, phi)
    # k = reduce_poly_mod(numerator * iden)

    # # Normalize and round the coefficients
    # k_coeffs = np.array(k.list())
    # k_coeffs = k_coeffs / int(a)
    # k_coeffs_rounded = np.round(k_coeffs).astype(int)

    # return R(list(k_coeffs_rounded[::-1]))

def key_generation(N, q):
    sigma_f = 1.17 * np.sqrt(q / (2 * N))
    
    i = 1
    while True:
        print(f"i: {i}")
        f = sample_polynomial(N, sigma_f)
        g = sample_polynomial(N, sigma_f)
        # print(f"f: {f}")
        # print(f"g: {g}")
        
        norm = gram_schmidt_norm(f, g)
        print(f"norm: {norm}")
        # pgcd, alpha, beta, rho_f, rho_g = pair_gcd(f, g)
        # print(f"pgcd: {pgcd}")
        valid, u, v, rho_f, rho_g = pair_gcd(f, g)

        if norm <= 1.17 * np.sqrt(q) and valid:
            print("norm and pgcd ok!")
            break
        i += 1

    print(f"u: {u}")
    print(f"v: {v}")
    print(f"rho_g: {rho_g}")

    F = -q * v * rho_g
    G = q * u * rho_f
    # print(f"F: {F}")
    # print(f"G: {G}")

    k = compute_k(f, g, F, G)



    return f, g

if __name__ == "__main__":
    # Generate key
    f, g = key_generation(N, q)