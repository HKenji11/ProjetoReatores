import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sb

# Estimativa de massa do catalisador e volume

densAparente =  1400 #kg/m³
L = 1.5 #m
d_tubo = 30e-2 #m

Vlamp = (0.0254/2)**2 * np.pi * L 
print(f"Volume de uma lâmpada {Vlamp * 1000} L")

#Volume de 1 reator com 1 lâmpada
Vreator = (d_tubo**2 * np.pi * L / 4 )  #Convertendo m³ para L
Vreator_verdadeiro = (Vreator - Vlamp) * 0.35

print(f"Volume de uma tubo {Vreator_verdadeiro * 1000} L")

Vcatalisador = (Vreator - Vlamp) * 0.65  #1 - 0.35 onde 0.35 é a fração de vazios
Mcatalisador = Vcatalisador * densAparente #kg


print(f"A massa de catalisador para uma lâmpada {Mcatalisador} kg")

# Declarando valores iniciais

## 2A -> C + D ## B -> C + E

#A -> H2O B -> Metanol C -> Hidrogênio D -> Oxigênio E-> Formaldeído

#condicoes de processo
P = 476.16 * 10 ** 3 #Pa
T = 150 + 273.15 #Kelvin

#inicialmente, vamos supor as vazões iniciais
fa0 = 2123 * 10        #mols/h
fb0 = 970 * 10         #mols/h

fc0 = fd0 = fe0 = 0     #mols/h
ft0 = fa0 + fb0         #mols/h
v0 = ft0 * 8.314 * T/P  #m³

ya0 = fa0/ft0
yb0 = fb0/ft0
yc0 = fc0/ft0
yd0 = fd0/ft0
ye0 = fe0/ft0

W = np.linspace(0, Mcatalisador*2, 10000)


# Declarando EDO's

def ODE(W, x):

    fa = x[0]
    fb = x[1]
    fc = x[2]
    fd = x[3]
    fe = x[4]
    ft = x[5]
    # ft = fa+fb+fc+fd+fe
    
    ya = fa/ft
    yb = fb/ft
    yc = fc/ft
    yd = fd/ft
    ye = fe/ft

    ra = (2/3) * ((1.7 * ya)  +  36.276*(yb/(0.038+ye*(1+ye/0.886))) )
    rb = (1/3) * ((1.7 * ya)  +  36.276*(yb/(0.038+ye*(1+ye/0.886))) )
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

parametros_iniciais = [fa0, fb0, fc0, fd0, fe0, ft0]
ans = sp.integrate.solve_ivp(ODE, t_span = (W[0], W[-1]), y0 = parametros_iniciais, t_eval = W)
ans

fb_list = ans.y[1,:]
fa_list = ans.y[0,:]
fc_list = ans.y[2,:]
xa = []
xb = []


for i in range(len(fb_list)):
    oi = (fb0 - fb_list[i])/fb0
    xb.append(oi)

for i in range(len(fa_list)):
    oi = (fa0 - fa_list[i])/fa0
    xa.append(oi)

print(fc_list)
print(xb)

kilos = fc_list[-1]*2

plt.plot(W,xb, color='green')
plt.plot(W,xa, color = "blue")
plt.show()