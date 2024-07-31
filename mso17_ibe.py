import numpy as np
from scipy.stats import norm as scipy_norm
from numpy.fft import fft, ifft

# Constants
N = 512
q = 2**30
sigma_1 = 0.84932180028801904272150283410288961971514109378435394286159953238339383120795466719298223538163406787061691601172910413284884326532697308797136114023
LDRMX = 2**31 - 1
log_2 = np.log(2)
omega = np.exp(2j * np.pi / N)
omega_1 = np.exp(-2j * np.pi / N)

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

def gram_schmidt_norm(f, g, N, q):
    # Norm of (g, -f)
    norm_1 = np.sqrt(np.sum(f**2 + g**2))
    print(f"norm_1: {norm_1}")

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
    print(f"norm_2: {norm_2}")

    return max(norm_1, norm_2)
    

def key_generation(N, q):
    sigma_f = 1.17 * np.sqrt(q / (2 * N))
    
    while True:
        f = sample_polynomial(N, sigma_f)
        g = sample_polynomial(N, sigma_f)
        # print(f"f: {f}")
        # print(f"g: {g}")
        
        norm = gram_schmidt_norm(f, g, N, q)


        if norm <= 1.17 * np.sqrt(q):
            break

    return f, g

if __name__ == "__main__":

    # Generate key
    f, g = key_generation(N, q)