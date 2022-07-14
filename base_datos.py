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

}
