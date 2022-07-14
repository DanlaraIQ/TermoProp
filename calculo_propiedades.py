# Archivo calcular las propiedades termodinámica de diferentes compuestos
import base_datos as bd
import math
import numpy as np
from scipy.integrate import quad

compuesto = "propano"


def calP_vap(compuesto, temperatura):
    '''
    En esta función se calcula a través del uso de la ecuación de Antoine la presión de vapor.
    Ecuación fue tomada de la base de datos de YAWS logP=A-(B/(T+C)). La presión está en milímetros de mercurio y la temperatura en Kelvin
    '''
    Apar = bd.propiedades[compuesto][0]
    Bpar = bd.propiedades[compuesto][1]
    Cpar = bd.propiedades[compuesto][2]
    P_vap = 10**(Apar - Bpar / (temperatura + Cpar - 273))
    P_vap = P_vap * 0.00133322  # convertir de milímetros de mercurio a bares
    return P_vap


def peng_robinson(compuesto, temperatura, presion):
    '''
    Esta función emplea la ecuación de Peng-Robinson para realizar el cálculo del factor de compresibilidad de líquido
    y vapor para el componente seleccionado además produce como salidas otras variables Termodinamicas que son útiles
    para diversos cálculos
    '''
    T_c = bd.propiedades[compuesto][6]
    P_c = bd.propiedades[compuesto][7]
    fac_ace = bd.propiedades[compuesto][11]
    T_r = temperatura / T_c
    P_r = presion / P_c
    m = 0.37464 + 1.5422 * fac_ace - 0.26992 * fac_ace**2
    alpha = (1 + m * (1 - math.sqrt(T_r)))**2
    A = 0.45724 * (P_r / T_r**2) * alpha
    B = 0.0778 * (P_r / T_r)

    polinomio = [1, -1 + B, A - 2 * B - 3 * B**2, -A * B + B**2 + B**3]
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
    return Z_vapor, Z_liquido, A, B, m, alpha


def cal_volEsp_tabla(compuesto, temperatura, presion, cte_R):
    '''
    Esta función calcula el volumen específico del líquido y vapor en metros cúbicos por kilogramo.
    '''
    peso_molecular = bd.propiedades[compuesto][5]
    Z_vapor, Z_liquido, A, B, m, alpha = peng_robinson(compuesto, temperatura, presion)
    if Z_liquido == None:
        volumen_especifico_vapor = (Z_vapor * cte_R * temperatura / presion) * (1000 / peso_molecular)
        volumen_especifico_liquido = "No hay líquido"
        return volumen_especifico_vapor, volumen_especifico_liquido
    else:
        volumen_especifico_liquido = (Z_liquido * cte_R * temperatura / presion) * (1000 / peso_molecular)
        volumen_especifico_vapor = (Z_vapor * cte_R * temperatura / presion) * (1000 / peso_molecular)
        return volumen_especifico_vapor, volumen_especifico_liquido


def cal_entalVaporizacion(compuesto, temperatura, cte_gases=8.314):
    '''
    Esta función se presenta la ecuación de Pitzer para el cálculo de la entalpía de evaporización
    Para una temperatura determinada
    '''
    T_c = bd.propiedades[compuesto][6]
    fac_ace = bd.propiedades[compuesto][11]
    peso_molecular = bd.propiedades[compuesto][5]
    T_r = temperatura / T_c
    cal_1 = 7.08 * (1 - T_r)**(0.354)
    cal_2 = 10.95 * (fac_ace) * (1 - T_r)**(0.456)
    entalpia_vaporizacion = ((cal_1 + cal_2) * (cte_gases * T_c)) / peso_molecular
    return entalpia_vaporizacion


def cal_base_entalVapor(compuesto, temperatura, presion):
    '''
    Esta función calcula la entalpía de vapor.
    '''
    peso_molecular = bd.propiedades[compuesto][5]
    T_c = bd.propiedades[compuesto][6]
    P_c = bd.propiedades[compuesto][7]
    fac_ace = bd.propiedades[compuesto][11]
    T_r = temperatura / T_c
    P_r = presion / P_c
    Z_vapor, Z_liquido, A, B, m, alpha = peng_robinson(compuesto, temperatura, presion)
    gamma = (0.37464 + 1.5422 * fac_ace - 0.26992 * fac_ace**2) * (math.sqrt(T_r / alpha))
    cal1 = Z_vapor - 1 - (A * (1 + gamma) / (math.sqrt(8) * B)) * \
        math.log((Z_vapor + (1 + math.sqrt(2)) * B) / (Z_vapor + (1 - math.sqrt(2)) * B))
    entalpia_vapor = cal1 * 8.314 * temperatura
    return entalpia_vapor


def cal_base_entroVapor(compuesto, temperatura, presion):
    '''
    Esta función calcula la entropia de vapor.
    '''
    peso_molecular = bd.propiedades[compuesto][5]
    T_c = bd.propiedades[compuesto][6]
    P_c = bd.propiedades[compuesto][7]
    fac_ace = bd.propiedades[compuesto][11]
    T_r = temperatura / T_c
    P_r = presion / P_c
    Z_vapor, Z_liquido, A, B, m, alpha = peng_robinson(compuesto, temperatura, presion)
    gamma = (0.37464 + 1.5422 * fac_ace - 0.26992 * fac_ace**2) * (math.sqrt(T_r / alpha))
    cal1 = math.log(Z_vapor - B) - ((A * gamma) / (math.log(8) * B)) * \
        math.log((Z_vapor + (1 + math.sqrt(2)) * B) / (Z_vapor + (1 - math.sqrt(2)) * B))
    entropia_vapor = cal1 * 8.314
    return entropia_vapor


def fun_Cp(T):
    Acp = bd.propiedades[compuesto][12]
    Bcp = bd.propiedades[compuesto][13]
    Ccp = bd.propiedades[compuesto][14]
    Dcp = bd.propiedades[compuesto][15]
    Ecp = bd.propiedades[compuesto][16]
    return Acp + Bcp * T + Ccp * T**2 + Dcp * T**3 + Ecp * T**4


def fun_Cp_entropia(T):
    Acp = bd.propiedades[compuesto][12]
    Bcp = bd.propiedades[compuesto][13]
    Ccp = bd.propiedades[compuesto][14]
    Dcp = bd.propiedades[compuesto][15]
    Ecp = bd.propiedades[compuesto][16]
    return (Acp + Bcp * T + Ccp * T**2 + Dcp * T**3 + Ecp * T**4) / T


def cal_entalpia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref):
    peso_molecular = bd.propiedades[compuesto][5]
    intCp = quad(fun_Cp, T_ref, temperatura)
    ental_1 = cal_base_entalVapor(compuesto, T_ref, P_ref)
    ental_2 = cal_base_entalVapor(compuesto, temperatura, presion)
    entalpia_referencia = cal_entalVaporizacion(compuesto, T_ref)
    entalpia_vapor = entalpia_referencia + (-ental_1 + intCp[0] + ental_2) / peso_molecular
    return entalpia_vapor


def cal_entropia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref):
    peso_molecular = bd.propiedades[compuesto][5]
    intCp_entro = quad(fun_Cp_entropia, T_ref, temperatura)
    intCp_entro = intCp_entro[0] / peso_molecular
    entropia_vaporizacion_ref = cal_entalVaporizacion(compuesto, T_ref) / temperatura
    entropia_vaporizacion = cal_entalVaporizacion(compuesto, temperatura) / temperatura
    entro_vapor_1 = cal_base_entroVapor(compuesto, T_ref, P_ref) / peso_molecular
    entro_vapor_2 = cal_base_entroVapor(compuesto, temperatura, presion) / peso_molecular
    entropia_vapor = entropia_vaporizacion_ref + (entro_vapor_2 - entro_vapor_1 + intCp_entro)
    return entropia_vapor


def cal_entropia_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref):
    peso_molecular = bd.propiedades[compuesto][5]
    entropia_vaporizacion = cal_entalVaporizacion(compuesto, temperatura) / temperatura
    entropia_vapor = cal_entropia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    entropia_liquido = entropia_vapor - entropia_vaporizacion
    return entropia_liquido


def cal_entropia_vaporizacion_tabla(compuesto, temperatura, presion, T_ref, P_ref):
    peso_molecular = bd.propiedades[compuesto][5]
    entropia_vaporizacion = cal_entalVaporizacion(compuesto, temperatura) / temperatura
    return entropia_vaporizacion


def cal_entalpia_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref):
    peso_molecular = bd.propiedades[compuesto][5]
    entalpia_vapor = cal_entalpia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    entalpia_referencia = cal_entalVaporizacion(compuesto, temperatura)
    entalpia_liquido = entalpia_vapor - entalpia_referencia
    return entalpia_liquido


def cal_enInterna_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref):
    peso_molecular = bd.propiedades[compuesto][5]
    Z_vapor, Z_liquido, A, B, m, alpha = peng_robinson(compuesto, temperatura, presion)
    entalpia_vapor = cal_entalpia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    energia_interna_vapor = entalpia_vapor - (temperatura * 8.314 * Z_vapor) / peso_molecular
    return energia_interna_vapor


def cal_enInterna_vaporizacion_tabla(compuesto, temperatura, presion, T_ref, P_ref):
    peso_molecular = bd.propiedades[compuesto][5]
    Z_vapor, Z_liquido, A, B, m, alpha = peng_robinson(compuesto, temperatura, presion)
    entalpia_vaporizacion = cal_entalVaporizacion(compuesto, temperatura)
    energia_interna_vaporizacion = entalpia_vaporizacion - (8.314 * temperatura * (Z_vapor - Z_liquido)) / peso_molecular
    return energia_interna_vaporizacion


def cal_enInterna_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref):
    energia_interna_vapor = cal_enInterna_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    energia_interna_vaporizacion = cal_enInterna_vaporizacion_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    energia_interna_liquido = energia_interna_vapor - energia_interna_vaporizacion
    return energia_interna_liquido
