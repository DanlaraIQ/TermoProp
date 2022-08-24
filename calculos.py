'''
    Este documento incluye las funciones necesarias para realizar los cálculos de las propiedades termodinámicas.
    La ecuación de Antoine es log10(P) = A + B/T + C*log10(T) + D*T + E*T^2
'''
import base_datos as bd

import numpy as np
import math


def cal_cp(T_actual):
    par = bd.cp_data[componente]
    return par[0] + par[1] * T_actual + par[2] * T_actual**2 + par[3]


def antoine(componente, T_actual):
    valores = bd.propiedades[componente]
    A = valores[1]
    B = valores[2]
    C = valores[3]
    presion_vapor = 10**(A - B / (T_actual + C))
    return presion_vapor


def van_der_Waals(componente, P_actual, T_actual):
    valores = bd.propiedades[componente]
    T_critica = valores[4]
    P_critica = valores[5]
    P_red = P_actual / P_critica
    T_red = T_actual / T_critica
    A = (27 / 64) * (P_red / T_red**2)
    B = (1 / 8) * (P_red / T_red)
    polinomio = [1, -1 - B, A, -A * B]
    sol_Z = np.roots(polinomio)
    Z = []
    for i in range(len(sol_Z)):
        if np.isreal(sol_Z[i]) == 1:
            Z.append(sol_Z[i].real)
    if len(Z) == 3:
        Z_vapor = max(Z)
        Z_liquido = min(Z)
    elif len(Z) == 2:
        Z_vapor = max(Z)
        Z_liquido = min(Z)
    elif len(Z) == 1:
        Z_vapor = max(Z)
    else:
        print("No se encontraron raíces positivas")
    m, alfa = 0, 0
    return Z_vapor, Z_liquido, A, B, m, alfa


def reclich_kwong(componente, P_actual, T_actual):
    valores = bd.propiedades[componente]
    T_critica = valores[4]
    P_critica = valores[5]
    P_red = P_actual / P_critica
    T_red = T_actual / T_critica
    A = 0.42748 * (P_red / T_red**2.5)
    B = 0.08664 * (P_red / T_red)
    polinomio = [1, -1, A - B - B**2, -A * B]
    sol_Z = np.roots(polinomio)
    Z = []
    for i in range(len(sol_Z)):
        if np.isreal(sol_Z[i]) == 1:
            Z.append(sol_Z[i].real)
    if len(Z) == 3:
        Z_vapor = max(Z)
        Z_liquido = min(Z)
    elif len(Z) == 2:
        Z_vapor = max(Z)
        Z_liquido = min(Z)
    elif len(Z) == 1:
        Z_vapor = max(Z)
        Z_liquido = None
    else:
        print("No se encontraron raíces positivas")
    m, alfa = 0, 0
    return Z_vapor, Z_liquido, A, B, m, alfa


def soave_reclich_kwong(componente, P_actual, T_actual):
    valores = bd.propiedades[componente]
    T_critica = valores[4]
    P_critica = valores[5]
    fac_ace = valores[8]
    P_red = P_actual / P_critica
    T_red = T_actual / T_critica
    m = 0.48 + 1.574 * fac_ace - 0.176 * fac_ace**2
    alfa = (1 + m * (1 - math.sqrt(T_red)))**2
    A = 0.42748 * (P_red / T_red**2) * alfa
    B = 0.08664 * (P_red / T_red)
    polinomio = [1, -1, A - B - B**2, -A * B]
    sol_Z = np.roots(polinomio)
    Z = []
    for i in range(len(sol_Z)):
        if np.isreal(sol_Z[i]) == 1:
            Z.append(sol_Z[i].real)
    if len(Z) == 3:
        Z_vapor = max(Z)
        Z_liquido = min(Z)
    elif len(Z) == 2:
        Z_vapor = max(Z)
        Z_liquido = min(Z)
    elif len(Z) == 1:
        Z_vapor = max(Z)
        Z_liquido = None
    else:
        print("No se encontraron raíces positivas")
    return Z_vapor, Z_liquido, A, B, m, alfa


def peng_robinson(componente, P_actual, T_actual):
    valores = bd.propiedades[componente]
    T_critica = valores[4]
    P_critica = valores[5]
    fac_ace = valores[8]
    P_red = P_actual / P_critica
    T_red = T_actual / T_critica
    m = 0.37464 + 1.5422 * fac_ace - 0.26992 * fac_ace**2
    alfa = (1 + m * (1 - math.sqrt(T_red)))**2
    A = 0.45724 * (P_red / T_red**2) * alfa
    B = 0.0778 * (P_red / T_red)
    polinomio = [1, -1 + B, A - 2 * B - 3 * B**2, -A * B + B**2 + B**3]
    sol_Z = np.roots(polinomio)
    print()
    Z = []
    for i in range(len(sol_Z)):
        if np.isreal(sol_Z[i]) == 1:
            Z.append(sol_Z[i].real)
    if len(Z) == 3:
        Z_vapor = max(Z)
        Z_liquido = min(Z)
    elif len(Z) == 2:
        Z_vapor = max(Z)
        Z_liquido = min(Z)
    elif len(Z) == 1:
        Z_vapor = max(Z)
        Z_liquido = None
    else:
        print("No se encontraron raíces positivas")
    return Z_vapor, Z_liquido, A, B, m, alfa


def vol_vapor(Z_vapor, cte_R, T_actual):
    peso_molecular = bd.propiedades[componente][0]
    vol_esp_vapor = (Z_vapor * cte_R * T_actual * 1000) / (P_actual * peso_molecular)
    return vol_esp_vapor


def vol_liquido(Z_vapor, cte_R, T_actual):
    peso_molecular = bd.propiedades[componente][0]
    vol_esp_liquido = (Z_liquido * cte_R * T_actual * 1000) / (P_actual * peso_molecular)
    return vol_esp_liquido


componente = "propano"

T_actual = 273 - 40  # K
P_actual = antoine(componente, T_actual)  # Pa
T_ref = 273
P_ref = antoine(componente, T_ref)  # Pa
cte_R = 8.3144626
peso_molecular = bd.propiedades[componente][0]

Z_vapor_ref, Z_liquido_ref, A_ref, B_ref, m_ref, alfa_ref = peng_robinson(componente, P_ref, T_ref)
print(Z_vapor_ref, Z_liquido_ref, A_ref, B_ref, m_ref, alfa_ref)
Z_vapor, Z_liquido, A, B, m, alfa = peng_robinson(componente, P_actual, T_actual)
volumen_vapor = vol_vapor(Z_vapor, cte_R, T_actual)
volumen_liquido = vol_vapor(Z_liquido, cte_R, T_actual)

B = 1274.89 * math.log(10)
C = 140.94
T_actual = 449.8
R = 8.314

h_vap = (B * R * T_actual**2) / (T_actual + C - 273.15)**2
print(h_vap)
