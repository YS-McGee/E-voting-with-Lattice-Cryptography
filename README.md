# Documentation

## Declaration
Lattice-IBE-master is referenced to tprest/Lattice-IBE [1], with modifications such as execution on Apple Silicon documented in the Issues [2], and additional functions called by main.py.

## Installation

SageMath Installation: https://doc.sagemath.org/html/en/tutorial/index.html

NTL.h and GMP.h on Apple Silicon:

```
$ brew install gmp
$ brew install ntl
```
## Execution

```$ sage -python mso17_ibe.py```

## Voting Parties

There are 3 parties invloved, Vd (voting device/voter), AT (voting authority), and Tk (trustees)

## Voting Procedures

### 1. Setup

AT: publishes public parametres -> N (IBE and Commitment lattice dimension), q (modulus), Ncand (number of candidate = 1), list of eligible voters, opening/closing time, and IDs for this election.
AT also runs KeyGeneration for the commitment scheme.

Each Tk runs KeyGeneration for identity-based encryption (IBE) and publishes Master Public Key (h).

### 2. Ballot Creation

Each Vd generates a vote of either 0 or 1 (agree/disagree), then use Shamir's secret sharing to creates a number of shares.




## Refereces

[1] https://github.com/tprest/Lattice-IBE/tree/master

[2] Compilation on Apple Silicon #2. https://github.com/tprest/Lattice-IBE/issues/2
