import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

# Componentes
# A:H2O 
# B:Metanol 
# C:Hidrogênio 
# D:Oxigênio 
# E:Formaldeído

# Declarando EDO para resolução dos balanços
def EDO(W, x):

    fa = x[0]
    fb = x[1]
    ft = x[5]
    
    ya = fa/ft
    yb = fb/ft

    ra = (2/3) * ((1.7 * ya)  +  36.276*(yb/(0.038+yb*(1+yb/0.886))) )
    rb = ra/2
    rc = ra + rb
    rd = ra/2
    re = rb
    dFAdW = -ra
    dFBdW = -rb
    dFCdW = rc
    dFDdW = rd
    dFEdW = re
    dFTdW = -ra - rb + rc + rd + re

    return [dFAdW, dFBdW, dFCdW, dFDdW, dFEdW, dFTdW]


# Dados iniciais do projeto
diamTubo = 20e-2        # Diâmetro interno do reator em metros
diamLamp = 0.026        # Diâmetro da lâmpada em metros
espessura = 0.001       # Espessura da camisa
comprimentoLamp = 1.5   # Comprimento da lâmpada em metros
espacamento = 0.001     # Espaçamento entre a lâmpada e a camisa
densAparente =  1400    # Densidade aparente do leito fixo em kg/m³
fracaoVazios = 0.35     # Fração de vazios


# Cálculo do volume da lâmpada
volumeLamp = ((diamLamp + espacamento + espessura)/2)**2 * np.pi * comprimentoLamp
print(f"Volume de uma lâmpada {volumeLamp * 1000} L")


# Cálculo do volume de um reator
Vreator = (diamTubo**2 * np.pi * comprimentoLamp / 4 )  #Convertendo m³ para L


# Cálculo do volume de um reator considerando a fração de vazios
volumeRealReator = (Vreator - volumeLamp) * fracaoVazios
print(f"Volume de um tubo {volumeRealReator * 1000} L")


# Cálculo do volume e massa do catalisador
volumeCat = (Vreator - volumeLamp) * (1 - fracaoVazios)  #1 - 0.35 onde 0.35 é a fração de vazios
massaCat = volumeCat * densAparente #kg
print(f"A massa de catalisador para uma lâmpada {massaCat} kg")


# Estimativa das vazões molares dos componentes Água e Metanol, em mol/h
fa0 = 700
fb0 = 846

fc0 = fd0 = fe0 = 0              
ft0 = fa0 + fb0 + fc0 + fd0 + fe0

ya0 = fa0/ft0
yb0 = fb0/ft0
yc0 = fc0/ft0
yd0 = fd0/ft0
ye0 = fe0/ft0

W = np.linspace(0, massaCat, 10000)

# Cálculo das vazões por método numérico - Scipy
vazoesIniciais = [fa0, fb0, fc0, fd0, fe0, ft0]
vazoesCalc = sp.integrate.solve_ivp(EDO, t_span = (W[0], W[-1]), y0 = vazoesIniciais, t_eval = W)

faList = vazoesCalc.y[0,:]
fbList = vazoesCalc.y[1,:]
fcList = vazoesCalc.y[2,:]
fdList = vazoesCalc.y[3,:]
feList = vazoesCalc.y[4,:]


# Declaração de listas para as conversões de água e metanol
xa = []
xb = []


# Cáculo das conversões de água e metanol
for i in range(len(faList)):
    conversao_calc = (fa0 - faList[i])/fa0
    xa.append(conversao_calc)

for i in range(len(fbList)):
    conversao_calc = (fb0 - fbList[i])/fb0
    xb.append(conversao_calc)

print("A conversão de H2O é ", xa[-1])
print("A conversão de MeOH é ", xb[-1])


# Produção de H2
taxa_mols = fcList[-1]
taxa_kg = taxa_mols * 2.016/1000
print("A produção de H2 é ", fcList[-1], "mols/h =", taxa_kg, "kg/h")


# Cálculo das composições finais de cada componente
vazoes_finais = [faList[-1], fbList[-1], fcList[-1], fdList[-1], feList[-1]]
FT_final = 0
for i in vazoes_finais:
    FT_final += i

composicoes = []
for j in vazoes_finais:
    composicoes.append(j/FT_final)

print("A composição final de H2O é ", composicoes[0])
print("A composição final de MeOH é ", composicoes[1])
print("A composição final de H2 é ", composicoes[2])
print("A composição final de O2 é ", composicoes[3])
print("A composição final de Formaldeído é ", composicoes[4])


# Plotagem dos resultados por gráficos

# Gráfico de taxa de conversão
plt.title('Taxas de Conversão em função da Massa do Catalisador', fontweight='bold')
plt.plot(W,xb, linestyle='solid', color='black', label='MeOH')
plt.plot(W,xa, linestyle='--', color = "blue", label='H2O')
plt.xlabel('Massa do Catalisador [kg]')
plt.ylabel('Taxa de Conversão')
plt.legend(facecolor='white', framealpha = 1)
plt.axhline(y=xb[-1], color='black', linestyle='dotted', alpha=0.2)
plt.text(W[-1] - 5, xb[-1], f'{xb[-1]:.3f}', color='black', fontsize=8)
plt.show()

# Gráfico de vazões molares
plt.title('Vazões Molares em função da Massa do Catalisador', fontweight='bold')
plt.plot(W,fbList, linestyle='-.', color='purple', label='MeOH')
plt.plot(W,faList, linestyle='--', color='blue', label='H2O')
plt.plot(W,fcList, linestyle='solid', color='black', label='H2')
plt.plot(W,fdList, linestyle='solid', color='green', label='O2')
plt.plot(W,feList, linestyle=':', color='black', label='H2CO')
plt.xlabel('Massa do Catalisador [kg]')
plt.ylabel('Vazão Molar [mol/h]')
plt.legend(facecolor='white', framealpha = 1)
plt.show()