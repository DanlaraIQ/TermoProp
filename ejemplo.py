import numpy as np
import math
# ln(P_vap) = A + B/(C + temp)
# PV = nRT


def ejemplo(P, V, mol, R):
    T = P * V / (mol * R)
    return T


def antoine(A, B, C, temp):
    P_vap = math.exp(A + B / (C + temp))
    return P_vap


def llama_funcion(A, B, C, temp, P, V, mol, R):
    r = antoine(A, B, C, temp) * ejemplo(P, V, mol, R)
    return r
