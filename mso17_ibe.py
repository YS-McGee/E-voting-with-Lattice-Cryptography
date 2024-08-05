from decimal import Decimal, getcontext
import numpy as np
from numpy.fft import fft, ifft
from scipy.stats import norm as scipy_norm
from sympy import symbols, ZZ, Poly
import sys
from sage.all import *
import gc


np.set_printoptions(suppress=True)

""" Constants """
# Set the precision to a high value to accommodate the given precision
getcontext().prec = 100
# Declare PiPrime with the given value
# PiPrime = Decimal('0.39894228040143267793994605993438186847585863116493465766592582967065792589930183850125233390730693643030255886263518268551099195455583724299621273062')
PiPrime = 0.39894228040143267793994605993438186847585863116493465766592582967065792589930183850125233390730693643030255886263518268551099195455583724299621273062
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

# def reduce_poly_mod(poly):
#     reduced = [(poly[i] - sum(poly[j] for j in range(N, len(poly)) if j - i == N)) for i in range(N)]
#     return reduced

def reduce_poly_mod(poly):
    if not poly:
        print("poly is empty")
        return [0] * N

    reduced = []
    # print(f"poly")
    for i in range(N):
        reduction_sum = sum(poly[j] for j in range(N, len(poly)) if j - i == N)
        # print(f"poly: {poly}")
        # print(f"i: {i}")
        # print(f"{poly[i]}")
        reduced.append(poly[i] - reduction_sum)
         
    return reduced


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
    # print(f"poly1: {poly1}")
    # print(f"poly2: {poly2}")
    # print(f"result_coeffs: {result_coeffs}")
    # print(f"P_result: {P_result}")
    # print(f"result_coeffs: {result_coeffs}")
    
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
    # print(f"rho_f: {rho_f}")
    
    # Return the inverse polynomial
    return inv_f * rho_f

def mpk_gen(f, g):
    # Define the polynomial ring over the finite field GF(q)
    R = PolynomialRing(GF(q), 'x')
    x = R.gen()

    # print(f"f: {f}")
    # print(f"g: {g}")

    # f = R(f)
    # g = R(g)

    f_inv = Inverse(f)
    # print(f"f_inv: {f_inv}")
    f_inv = sage_to_numpy(f_inv).tolist()
    # print(f"f_inv: {f_inv}")
    # print(f"f: {f}")
    # print(f"g: {g}")
    # print(f"R(g): {R(g)}")
    # print(f"f_inv: {f_inv}")
    h = [a%q for a in multiply_large_polynomials(g, f_inv)]

    # print(f"h: {h}")
    return h

def anticirculant_matrix(f):
    N0 = len(f)
    A = np.zeros((N0, N0), dtype=int)
    
    for i in range(N0):
        # Fill in the first part directly from the list f
        A[i, i:] = f[:N0-i]
        # Fill in the second part with negative values adjusted for the wrap-around
        A[i, :i] = -np.array(f[N0-i:])
    
    return A

def b_matrix(f, g, F, G):
    neg_f = [-x for x in f]
    neg_F = [-x for x in F]

    # print(f"f: {f}")
    # print(f"g: {g}")
    # print(f"neg_f: {neg_f}")
    # print(f"neg_F: {neg_F}")

    M = np.zeros((2 * N, 2 * N), dtype=int)

    # Fill top-left submatrix with A(g)
    M[:N, :N] = anticirculant_matrix(g)
    
    # Fill top-right submatrix with -A(f)
    M[:N, N:] = anticirculant_matrix(neg_f)
    
    # Fill bottom-left submatrix with A(G)
    M[N:, :N] = anticirculant_matrix(G)
    
    # Fill bottom-right submatrix with -A(F)
    M[N:, N:] = anticirculant_matrix(neg_F)

    # print(M)
    return M

def key_generation():
    sigma_f = 1.17 * np.sqrt(q / (2 * N))
    
    # i = 1
    while True:
        # print(f"i: {i}")
        f = sample_polynomial(N, sigma_f)
        g = sample_polynomial(N, sigma_f)        
        norm = gram_schmidt_norm(f, g)
        valid, u, v, sage_rho_f, sage_rho_g = pair_gcd(f, g)

        if norm <= 1.17 * np.sqrt(q) and valid:
            # print(f"norm: {norm}")
            print("norm and gcd ok!")
            break
        # i += 1

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
    if q != q_test:
        print("FAILED f*G - g*F != q")
        sys.exit()

    mpk = mpk_gen(f, g)
    # print(f"f: {f}")
    # print(f"-f: {[-x for x in f]}")
    msk = b_matrix(f, g, F, G)

    # print(f"MPK: {MPK}")
    
    return mpk, msk, f, g

def fast_mgs(B):
    """
    Perform the Fast Modified Gram-Schmidt (MGS) orthogonalization on matrix B.
    """
    B = np.array(B, dtype=float)
    Bst = np.zeros_like(B)

    # Initial setup
    Bst[0, :] = B[0, :]
    v = np.zeros(2 * N)
    v[:N-1] = Bst[0, 1:N]
    v[N-1] = -Bst[0, 0]
    v[N:2*N-1] = Bst[0, N+1:2*N]
    v[2*N-1] = -Bst[0, N]
    v1 = v.copy()
    C_k = np.dot(Bst[0], v)
    D_k = np.dot(v, v)

    # Orthogonalize first half
    for k in range(1, N):
        aux = C_k / D_k
        Bst[k, 0] = -Bst[k-1, N-1] + aux * v[N-1]
        Bst[k, N] = -Bst[k-1, 2*N-1] + aux * v[2*N-1]
        for j in range(1, N):
            Bst[k, j] = Bst[k-1, j-1] - aux * v[j-1]
            Bst[k, j+N] = Bst[k-1, j+N-1] - aux * v[j+N-1]
        v -= aux * Bst[k-1]
        C_k = np.dot(Bst[k], v1)
        D_k -= (C_k * C_k) / D_k

    # Normalize and orthogonalize second half
    D_k = np.dot(Bst[N-1], Bst[N-1])
    for j in range(N):
        Bst[N, N+j] = Bst[N-1, N-1-j] * q / D_k
        Bst[N, j] = -Bst[N-1, 2*N-1-j] * q / D_k
    v[:N-1] = Bst[N, 1:N]
    v[N-1] = -Bst[N, 0]
    v[N:2*N-1] = Bst[N, N+1:2*N]
    v[2*N-1] = -Bst[N, N]
    v1 = v.copy()
    C_k = np.dot(Bst[N], v1)
    D_k = np.dot(Bst[N], Bst[N])

    for k in range(N + 1, 2 * N):
        aux = C_k / D_k
        Bst[k, 0] = -Bst[k-1, N-1] + aux * v[N-1]
        Bst[k, N] = -Bst[k-1, 2*N-1] + aux * v[2*N-1]
        for j in range(1, N):
            Bst[k, j] = Bst[k-1, j-1] - aux * v[j-1]
            Bst[k, j+N] = Bst[k-1, j+N-1] - aux * v[j+N-1]
        v -= aux * Bst[k-1]
        C_k = np.dot(Bst[k], v1)
        D_k -= (C_k * C_k) / D_k

    return Bst

def gpv(c, MSKD):
    # Beware to use copy(), otherwise it will affect c when using ci
    ci = c.copy()

    for i in range(2*N - 1, -1, -1):
        aux = MSKD['GS_Norm'][i]
        cip = np.dot(ci, MSKD['Bstar'][i]) / (aux**2)
        sip = MSKD['Sigma'] / aux
        zi = sample4(cip, sip*PiPrime)

        for j in range(2 * N):
            ci[j] -= zi * MSKD['B'][i][j]

    sk = [0] * (2 * N)
    for j in range(0, 2 * N):
        sk[j] = c[j] - ci[j]

    return sk

def ibe_extract(id, MSKD):
    # Initialize c with zeros
    c = np.zeros(2 * N)

    # Conversion and assignment
    for i in range(N):
        # c[i] = float(id[i])  # conv<double>(id[i]) in C++
        c[i] = id[i]
        c[i + N] = 0

    # Print the result
    # print("c:", c)

    sk = gpv(c, MSKD)

    sk[:N] = [c[i] - sk[i] for i in range(N)]
    sk[N:] = [-sk[i + N] for i in range(N)]

    # print(f"sk: {sk}")
    SK_id = [[0] * N for _ in range(2)]  # Initialize SK_id as a 2D list with zeros
    SK_id[0] = [sk[i] for i in range(N)]
    SK_id[1] = [sk[i + N] for i in range(N)]
    # print(f"SK_id: {SK_id}")

    return SK_id

def ibe_verify_key(SK_id, id, MSKD):
    f = MSKD['Prk'][0]
    g = MSKD['Prk'][1]

    t = id.copy()

    alpha = [SK_id[0][i] - t[i] for i in range(N)]
    bravo = multiply_large_polynomials(alpha, f)
    charlie = multiply_large_polynomials(g, SK_id[1])
    delta = [(bravo[i] + charlie[i])%q for i in range(N)]
    # print(f"all(x==0 for x in delta): {all(x==0 for x in delta)}")
    # print(f"delta: {delta}")
    # delta[1] = 3

    # Verify if it's a zero list
    return all(x==0 for x in delta)

def extract_test(id, MSKD):
    SK_id = ibe_extract(id, MSKD)

    if not ibe_verify_key(SK_id, id, MSKD):
        print("[FAIL] --- Key Verification Failed")
        sys.exit()
    # print("================= Key Extraction Test Successful! =================")

    return SK_id

def ibe_encrypt(message, id, mpk):
    e1 = (np.random.randint(0, 3, N) - 1).tolist()
    e2 = (np.random.randint(0, 3, N) - 1).tolist()
    romeo = (np.random.randint(0, 3, N) - 1).tolist()

    product_1 = multiply_large_polynomials(romeo, mpk)
    # print(f"romeo: {romeo}")
    # print(f"mpk: {mpk}")
    # print(f"product_1: {product_1}")
    uniform = [((a + b+ (q // 2)))%q-(q//2) for a, b in zip(product_1, e1)]
    # print(f"uniform: {uniform}")
    product_2 = multiply_large_polynomials(romeo, id)
    victor = [((a+b+(q//2)*c+(q//2))%q)-(q//2) for a, b, c in zip(product_2, e2, message)]

    # uniform_2 = [(a+e)%q for a, e in zip(product_1, e1)]
    # victor_2 = [(a+e+(q//2)*m)%q for a, e, m in zip(product_2, e2, message)]
    # victor_2 = [2**20 * (x // 2**20) for x in victor_2]
    # # print(f"uniform_2: {uniform_2}")
    # # print(f"victor_2: {victor_2}")


    # product_3 = multiply_large_polynomials(uniform_2, sk[1])
    # whiskey = [(a-b)//(q//2) for a, b in zip(victor_2, product_3)]
    # # print(f"whiskey: {whiskey}")

    cipher = [uniform, victor]
    
    return cipher

def ibe_decrypt(cipher, SK_id):
    # print(f"SK_id[1]: {SK_id}")
    message = multiply_large_polynomials(cipher[0], SK_id)

    # print(f"message: {message}")

    for i in range(N):
        message[i] = (cipher[1][i] - message[i]) % q
        message[i] = (message[i] + (q >> 2)) // (q >> 1)
        # print(f"message[i]: {message[i]}")
        message[i] %= 2
    
    return message

def encrypt_test(j, mpk, MSKD):
    id = [np.random.randint(0, q - 1) for _ in range(N)]

    SK_id = ibe_extract(id, MSKD)
    # print(f"SK)id: {SK_id}")
    # extract_test(id, MSKD)

    message = np.random.randint(0, 2, N).tolist()
    # print(f"message: {message}")

    cipher = ibe_encrypt(message, id, mpk)
    decrypted = ibe_decrypt(cipher, SK_id[1])

    # print(f"cipher: {cipher}")
    # print(f"decrypted: {decrypted}")
    # print(f"j: {j}")
    if message != decrypted:
        print(f"---- Decryption Failed at j:{j} ----")
        return False
    else:
        print(f"{j}: Encryption/Decryption Successful")
        return True
        # print(f"id: {id}")
        # print(f"mpk: {mpk}")
        # print(f"SK_id: {SK_id}")
        # print(f"message: {message}")
        # print(f"cipher: {cipher}")
        # print(f"decrypted: {decrypted}")

    # message = []
    # cipher = []
    # decrypted = []
    # print(f"message: {message}")
    # print(f"decrypted: {decrypted}")

        # sys.exit()
    # else:
    #     print("SUCCESSFULLY")



if __name__ == "__main__":
    # Generate key
    mpk, msk, f, g = key_generation()

    # print("================= MPK & MSK Generated =================")
    # print(f"MPK: {MPK}")
    # print(f"MSK: {MSK}")
    # print(msk)


    bstar = fast_mgs(msk)
    gs_norm = [0] * 2*N
    for i in range(0, 2*N):
        gs_norm[i] = np.linalg.norm(bstar[i, :])
    
    sigma = 2 * gs_norm[0]
    MSKD = {
        'Prk': [f, g],
        'B': msk,
        'Bstar':bstar,
        'GS_Norm':gs_norm,
        'Sigma':sigma
    }
    # id = [np.random.randint(0, q - 1) for _ in range(N)]
    # SK_id = extract_test(id, MSKD)
    # message = np.random.randint(0, 2, N).tolist()
    # cipher = ibe_encrypt(message, id, mpk)
    # decrypted = ibe_decrypt(cipher, SK_id[1])

    # if message != decrypted:
    #     print(f"---- Decryption Failed ----")

    for j in range(0, 1):
        # print(f"j: {j}")
        id = [np.random.randint(0, q - 1) for _ in range(N)]
        SK_id = extract_test(id, MSKD)
        valid = encrypt_test(j, mpk, MSKD)

        # if valid:
        #     for i in range(0, 30):
        #         message = np.random.randint(0, 2, N).tolist()
        #         cipher = ibe_encrypt(message, id, mpk)
        #         decrypted = ibe_decrypt(cipher, SK_id[1])

        #         if message != decrypted:
        #             print(f"---- Decryption Failed at i:{i} ----")
        del SK_id, id
        gc.collect()

    del MSKD, mpk, f, g, bstar, gs_norm, sigma
    gc.collect()

    