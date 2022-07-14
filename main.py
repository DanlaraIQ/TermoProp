import base_datos as bd
import calculo_propiedades as calP
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
'''
Este código calcula las siguientes propiedades termodinámicas:
- Volúmen específico de líquido
- Volúmen específico del vapor
- Energía interna de líquido, vapor y cambio de fase
- Entropía de líquido, vapor y cambio de fase
- Entalpía de líquido, vapor y cambio de fase
'''


#- Escriba el nombre del compuesto
compuesto = "propano"

#- En la variable "tipo" escribir True si se quieren obtener la propiedades termodinámicas en un intervalo de temperaturas
#- En la variable "tipo" escribir False si se quieren obtener la propiedades termodinámicas para una sola temperatura

tipo = True

# -- No se recomienda manipular estos parámetros. Son los datos para calcular el estado de referencia.
T_ref = 273 + bd.propiedades[compuesto][18]
cte_R = 8.314e-5
P_ref = calP.calP_vap(compuesto, T_ref)
propTer = ["Témperatura (K)", "Presión de saturación (bar)", "Volúmen de líquido (m^3/Kg)",
           "Volúmen de vapor (m^3/Kg)", "Energía interna de líquido (kJ/kg)", "Energía interna de vaporización (kJ/kg)",
           "Energía interna de vapor (kJ/kg)", "Entalpía de líquido (kJ/kg)",  "Entalpía de vaporización (kJ/kg)", "Entalpía de vapor (kJ/kg)",
           "Entropía de líquido (kJ/kg-K)",  "Entropía de vaporización (kJ/kg-K)", "Entropía de vapor (kJ/kg-K)"]
# ---------------------------------------------------------------------------------------------------

#-- Si se quieren las propiedades para un solo valor, escriba el valor a continuación en Kelvin:
valor_temperatura = 273 + 20  # Para una sola temperatura
#--------------------------
temperatura_inicial = 273 - 90
temperatura_final = 273 + 90
numero_valores_intermedios = 50


#-- Si se quieren las propiedades termodinámicas para un rango de temperaturas

if tipo == False:
    temperatura = valor_temperatura
    presion = calP.calP_vap(compuesto, temperatura)
    volumen_especifico_vapor, volumen_especifico_liquido = calP.cal_volEsp_tabla(compuesto, temperatura, presion, cte_R)
    entalpia_vaporizacion = calP.cal_entalVaporizacion(compuesto, temperatura)
    entalpia_vapor = calP.cal_entalpia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    entalpia_liquido = calP.cal_entalpia_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    entropia_vapor = calP.cal_entropia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    entropia_vaporicacion = calP.cal_entropia_vaporizacion_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    entropia_liquido = calP.cal_entropia_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    energia_interna_vapor = calP.cal_enInterna_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    energia_interna_vaporizacion = calP.cal_enInterna_vaporizacion_tabla(compuesto, temperatura, presion, T_ref, P_ref)
    energia_interna_liquido = calP.cal_enInterna_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref)

    print("Para una temperatura de: ", temperatura, " K")
    print("La presión de saturación es: {:.4f}".format(presion), " bar")
    print("El volúmen específico de líquido es: {:.4f}".format(volumen_especifico_liquido), " m^3/kg")
    print("El volúmen específico de vapor es: {:.4f}".format(volumen_especifico_vapor), " m^3/kg")
    print("La energía interna del líquido es: {:.4f}".format(energia_interna_liquido), " kJ/kg")
    print("La energía interna del vaporización es: {:.4f}".format(energia_interna_vaporizacion), " kJ/kg")
    print("La energía interna del vapor es: {:.4f}".format(energia_interna_vapor), " kJ/kg")
    print("La entalpía del líquido es: {:.4f}".format(entalpia_liquido), " kJ/kg")
    print("La entalpía de vaporización es: {:.4f}".format(entalpia_vaporizacion), " kJ/kg")
    print("La entalpía del vapor es: {:.4f}".format(entalpia_vapor), " kJ/kg")
    print("La entropia del líquido es: {:.4f}".format(entropia_liquido), " kJ/kg-K")
    print("La entropia de evaporación es: {:.4f}".format(entropia_vaporicacion), " kJ/kg-K")
    print("La entropia del vapor es: {:.4f}".format(entropia_vapor), " kJ/kg-K")
elif tipo == True:
    val_temp = [273 - 90, 200, 230, 260, 290, 320, 350, 273 + 90] # np.linspace(temperatura_inicial, temperatura_final, numero_valores_intermedios)
    tabla_propiedades = np.zeros([len(val_temp), len(propTer)])
    # tabla_propiedades[0, :] = np.asarray(propTer).copy()
    for i in range(len(val_temp)):
        temperatura = val_temp[i]
        presion = calP.calP_vap(compuesto, temperatura)
        volumen_especifico_vapor, volumen_especifico_liquido = calP.cal_volEsp_tabla(compuesto, temperatura, presion, cte_R)
        entalpia_vaporizacion = calP.cal_entalVaporizacion(compuesto, temperatura)
        entalpia_vapor = calP.cal_entalpia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
        entalpia_liquido = calP.cal_entalpia_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref)
        entropia_vapor = calP.cal_entropia_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
        entropia_vaporicacion = calP.cal_entropia_vaporizacion_tabla(compuesto, temperatura, presion, T_ref, P_ref)
        entropia_liquido = calP.cal_entropia_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref)
        energia_interna_vapor = calP.cal_enInterna_vapor_tabla(compuesto, temperatura, presion, T_ref, P_ref)
        energia_interna_vaporizacion = calP.cal_enInterna_vaporizacion_tabla(compuesto, temperatura, presion, T_ref, P_ref)
        energia_interna_liquido = calP.cal_enInterna_liquido_tabla(compuesto, temperatura, presion, T_ref, P_ref)

        propiedades_temp = [temperatura, presion, volumen_especifico_liquido, volumen_especifico_vapor, energia_interna_liquido, energia_interna_vaporizacion,
                            energia_interna_vapor, entalpia_liquido, entalpia_vaporizacion, entalpia_vapor, entropia_liquido, entropia_vaporicacion, entropia_vapor]
        tabla_propiedades[i, :] = propiedades_temp.copy()

    tabla = pd.DataFrame(tabla_propiedades)
    tabla.columns = propTer
    tabla.to_excel("propiedades_" + compuesto + ".xlsx", index=False)

    # Gráfica de temperatura-presion
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 1])
    plt.xlabel("Temperatura, [K]")
    plt.ylabel("Presión, [bar]")
    plt.savefig("P_T" + compuesto + ".pdf")
    plt.clf()

    # Gráfica de temperatura-volumen de vapor
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 2])
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 3])
    plt.xlabel("Temperatura, [K]")
    plt.ylabel("Volumen específico, [m^3/kg]")
    plt.legend(["Líquido", "Vapor"])
    plt.savefig("V_T" + compuesto + ".pdf")
    plt.clf()

    # Gráfica de temperatura-energía interna
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 4])
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 6])
    plt.xlabel("Temperatura, [K]")
    plt.ylabel("Energía interna, [kJ/kg]")
    plt.legend(["Líquido", "Vapor"])
    plt.savefig("U_T" + compuesto + ".pdf")
    plt.clf()

    # Gráfica de temperatura-entalpía
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 7])
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 9])
    plt.xlabel("Temperatura, [K]")
    plt.ylabel("Entalpía, [kJ/kg]")
    plt.legend(["Líquido", "Vapor"])
    plt.savefig("H_T" + compuesto + ".pdf")
    plt.clf()

    # Gráfica de temperatura-entropía
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 10])
    plt.plot(tabla_propiedades[:, 0], tabla_propiedades[:, 12])
    plt.xlabel("Temperatura, [K]")
    plt.ylabel("Entropía, [kJ/kg-K]")
    plt.legend(["Líquido", "Vapor"])
    plt.savefig("S_T" + compuesto + ".pdf")
