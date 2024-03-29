'''
    Este documento contiene las propiedades termodinámicas necesarias para los cálculos de las propiedades.
    La ecuación de Antoine es log10(P) = A - B/(T + C - 273.15)
    Presión en mmHg
    B liq = a(1 - T/Tc)^m
    Temperatura en K
'''
# PM = 0
# A = 1
# B = 2
# C = 3
# Tc = 4
# Pc = 5
# Vc = 6
# den = 7
# Zc = 8
# fac ace = 9
propiedades = {
    # "nombre_componente": [A, B, C, Tmin Tmax MW Tc(K)  Pc(bar)  Vc(cm3/mol)  ρ(g/cm3)  Zc  ω [(A    B   C   D   E)-Cp] T_ref]
    "propano": [7.01887, 889.864, 257.084, -187.69, 96.67, 44.097, 369.82, 42.49, 202.9, 0.2174, 0.28, 0.152, 28.277, 0.11600, 1.9597E-04, -2.327E-07, 6.8669E-11, 9.9502E-04, -40],
    "isopentano": [7.03015, 1140.45, 247.012, -159.9, 187.28, 72.15, 460.43, 33.81, 305.8, 0.2359, 0.27, 0.228, -0.881, 0.47498, -2.4797E-04, 6.751E-08, -8.5343E-12, -40],
    "abietic acid": [7.06195, 1850.91, 66.1306, 173.5, 558.85,  302.457, 832, 16.8, 930, 0.3252,  0.226,   1.129, -260.14, 2.72960, -2.3498E-03, 1.024E-06,   -1.8203E-10, -40],
    "acenaphthene": [7.16063, 2033.35, 197.714, 93.41, 530, 154.211, 803.15,  31,  553, 0.2789,  0.257,   0.381,   -61.063, 0.96388, -7.2100E-04, 2.623E-07,   -3.6857E-11, -40],
    "acetal": [7.28439, 1445.63, 224.685, -100,    267.85,  118.176, 541, 29.8,    402, 0.294,   0.266,   0.432,   31.834, 0.34539, 5.6817E-04,  -1.027E-06,  4.4628E-10,  -40],
    "acetaldehyde": [7.25504, 1110.4,  233.451, -123,    187.85,  44.053,  461, 55.5,   157, 0.2806,  0.227,   0.317,   34.14,   0.04002, 1.5634E-04,  -1.644E-07,  4.7248E-11,  -40],
    "acetamide": [7.40993, 1808.24, 178.098, 81,  487.85,  59.068,  761, 66, 215, 0.2747,  0.224,   0.189,   17.748,  0.13627, 1.0668E-04,  -1.865E-07,  6.2842E-11,  -40]
}
