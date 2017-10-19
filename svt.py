import numpy as np
import operator
import matplotlib
import matplotlib.pyplot as plt
from numpy import linalg as la
import scipy.sparse as sps


# % svt : Singular Value Thresholding Algorithm
# % Description : Recover a low-rank matrix M from a subset of sampled entries
#
# % omega : sampled set
# % P_omegaM : sampled entries
# % delta : step size
# % epsilon : tolerance
# % tau : parameter (threshold)
# % el : increment
# % kMAX : maximum iteration count


def stripCols(matrix, n):
    m = []
    for row in matrix:
        m.append(list(row[:n]))
    return m

def hermitian(matrix):
    """
    creates the hermitian matrix of given matrix (list or numpy.ndarray)
    :param matrix: given matrix
    :return: hermitian matrix
    """

    # Row index
    i = 0
    for row in matrix:
        # Column index
        j = 0
        for e in row:
            matrix[i][j] = np.conjugate(e)
            j+=1
        i+=1
    return matrix


def svt(n,omega, P_omegaM, delta, epsilon, tau, el, kMAX):


    # % Initialize r & s integers for loop:
    # % *  Notated as vectors : Each element's index
    # % *  is the correlated to the k iteration for
    # % *  main for loop
    r = np.zeros(kMAX)
    s = np.zeros(kMAX)

    # % Initialize X because used in Y definition
    X = np.zeros(n[0], n[1])

    # Initialize X because used in Y definition
    RAN = np.zeros(kMAX)
    Rel = np.zeros(kMAX)
    Res = np.zeros(kMAX)
    best = np.zeros(kMAX)

    # % Norm of P_omega(M)
    normProjM = la.norm(P_omegaM)

    # Define Y(0)
    k_0 = np.ceil(tau / (delta * normProjM))
    Y = k_0*delta*P_omegaM

    # % r integer for first iteration
    R = 0

    for k in range(kMAX-1):
        iter = k
        s[k] = R + 1

        while True:
            # %always nonincreading    % When used P_omega(M), r only = 0))
            [U,S,V] = la.svd(Y,full_matrices=False)
            if S[s[k]][s[k]] <= tau:
                break
            s[k] = s[k] + 1

        # % Set r to the largest index of the singular value
        # % that is greater than tau
        sigma = np.array(np.diag(S))

        # % Set r to the largest index of the singular value that is greater than tau
        for SS in range(s[k]):
            if sigma[SS] > tau:
                r[k] = SS

        # r(k-1) for kth iteration
        R = r[k]

        # Resizing SVD pieces
        U = stripCols(U, r[k])
        S = stripCols(S[:r[k]], r[k])
        V = stripCols(V, r[k])















def main():
    pass












