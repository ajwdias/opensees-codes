import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import pandas as pd 
import math as math

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Geometria do pilar
B = 0.55 #m
H = 0.55 #m
L = 1.65 #m
A = B*H   # área da seção transversal
I = (B*H**3)/12 # inércia

#uniaxialMaterial('Steel01', matTag, Fy, E0, b, a1, a2, a3, a4)
matTag_aco = 1
Fy = 511e6 #N/m2 resistência ao escoamento
E0 = 200e9 #N/m2 mod de elasticidade
b = 0.01
ops.uniaxialMaterial('Steel01', matTag_aco, Fy, E0, b)

#NAO CONFINADO
#uniaxialMaterial('Concrete01', matTag, fpc, epsc0, fpcu, epsU)
matTag_NC = 2
fpc = 32e6 #N/m2
epsc0 = 0.00208 #E0 = 22.9 N/m2
fpcu = 0.2 * fpc #aproximado 
epsU = 0.0045 #aproximado do grafico
ops.uniaxialMaterial('Concrete01', matTag_NC, fpc, epsc0, fpcu, epsU)


#Variáveis para calcular concreto confinado 1
epsc = epsc0
Ec = 2 * fpc / epsc0 #tg do grafico do concrete01
cobrimento = 0.04 #m
barras_camadas_1 = 5
diam_transversal_1 = 0.012 
diam_longitudinal_1 = 0.020
npernas_yy = npernas_xx = 4
espacamento_s_1 = 0.11
b_beam = 0.55
h_beam = 0.55
fy_transv_1 = 325 #MPa
fc_inicial = 32e6 
ConcretoC_Tag_1 = 3



# Calculation of Unconfined Concrete Strength and Strain

area_longitudinal = diam_longitudinal_1**2 * math.pi / 4
area_total_separada = n_barras_i * area_longitudinal for n_barras_i in barras_camadas_1
Ast = sum(area_total_separada)

# area_longitudinal = diam_longitudinal_1**2 * math.pi / 4
# #area_longitudinal = [diametro**2 * math.pi / 4 for diametro in diam_longitudinal_1]
# area_total_separada = [n_barras_i * area_i for n_barras_i, area_i in zip(barras_camadas_1, area_longitudinal)]
# Ast = sum(area_total_separada)
    # area_longitudinal = 3.1415/4*self.diam_longitudinal**2
    # n_barras = sum(self.barras_camadas)
    # Ast = n_barras * area_longitudinal
    
n_camadas = len(barras_camadas_1)
      
dlinha = cobrimento + diam_transversal_1 + diam_longitudinal_1/2
    # dlinha = self.cobrimento + self.diam_transversal + self.diam_longitudinal/2
slinha = espacamento_s_1 - diam_transversal_1
    
hx = b_beam - 2*cobrimento - diam_transversal_1   # dc
hy = h_beam - 2*cobrimento - diam_transversal_1   # bc
    
wlinha_x = ((hx - diam_transversal_1 - 2*diam_longitudinal_1)-(npernas_yy-2)*diam_longitudinal_1) / (npernas_yy-1)
wlinha_y = ((hy - diam_transversal_1 - 2*diam_longitudinal_1)-(npernas_xx-2)*diam_longitudinal_1) / (npernas_xx-1)

num_wlinha_x = 2 * (npernas_yy - 1)  
num_wlinha_y = 2 * (npernas_xx - 1) 
    
rho_cc = Ast/(hx * hy)  # ratio of area of longitudinal reinforcement to area of core of section 
Ae = (hx * hy - num_wlinha_x*(wlinha_x**2)/6 - num_wlinha_y*(wlinha_y**2)/6)*(1-0.5*slinha/hx)*(1-0.5*slinha/hy)   # Area of effectively confined concrete core
Acc = (hx * hy) * (1 - rho_cc) # Area of core concrete without longitudinal reinforcement (le resta área al núcle de concreto confinado hx*hy ó dc*bc)
ke = Ae / Acc    # Coefinement effectiveness coefficient
    
rho_x = (npernas_xx * math.pi * diam_transversal_1**2 / 4) / (espacamento_s_1*hy)   # Total area of transverse reinforcement parallel to the x axis
rho_y = (npernas_yy * math.pi * diam_transversal_1**2 / 4) / (espacamento_s_1*hx)   # Total area of transverse reinforcement parallel to the y axis
    
flx = ke * rho_x * fy_transv_1 * (10**6) # Lateral confining stress on the concrete in the direction x
fly = ke * rho_y * fy_transv_1 * (10**6) # Lateral confining stress on the concrete in the direction y
q = np.minimum(flx,fly)/np.maximum(flx,fly) 
    
A = 6.8886 - (0.6069 + 17.275*q)*np.exp(-4.989*q)
B = ((4.5)/(5/A*(0.9849-0.6306*np.exp(-3.8939*q))-0.1))-5
xlinha = (flx + fly) / (2*(fc_inicial))
k1 = A*(0.1+0.9/(1+B*xlinha))
k2 = 5*k1  # para resistencia normal da armadura transversal (5k1)
# k2 = 3*k1  # para elevada resistencia da armadura transversal (3k1)
fator_fc_confinado = 1+k1*xlinha
fator_ec_confinado = 1+k2*xlinha
fc_confinado = fc_inicial*fator_fc_confinado          #Resistência máxima do concreto confinado à compressão
ec_confinado = epsc*fator_ec_confinado         #Deformação de pico do concreto confinado à compressão 
if fc_confinado > fc_inicial:
    xn = 30                                                      #ok
    n = Ec*ec_confinado/(fc_confinado)                           #ok
    r = n/(n-1)                                                  #ok
else:
    fc_confinado = fc_inicial
    ec_confinado = epsc
xn = 2.3
    


epsc_u_NC2 = 0.008  #0.008 La defomración unitaria de ruptura del hormigón CONFINADO no es dado en el paper 34. Se asume 0.008. 

    # ops.uniaxialMaterial('Concrete01', matTag, fpc, epsc0, fpcu, epsU)
ops.uniaxialMaterial('Concrete01', ConcretoC_Tag_1, -fc_confinado, -ec_confinado, -fc_confinado*0.20, -epsc_u_NC2)  

#CONFINADO 1 - 55cm
matTag_C = ConcretoC_Tag_1
fpc_c = -fc_confinado
epsc0_c = -ec_confinado
fpcu_c = -fpc_u_NC
epsU_c = -epsc_u_NC2
ops.uniaxialMaterial('Concrete01', matTag_C, fpc_c, epsc0_c, fpcu_c, epsU_c)



#Variáveis para calcular concreto confinado 1
cobrimento = 0.04 #m
barras_camadas_2 = 5
diam_transversal_2 = 0.012 
diam_longitudinal_2 = 0.020
npernas_yy = npernas_xx = 4
espacamento_s_2 = 0.22
b_beam = 0.55
h_beam = 0.55
fy_transv_2 = 325 #MPa
fc_inicial = 32e6
ConcretoC_Tag_2 = 3



# Calculation of Unconfined Concrete Strength and Strain
area_longitudinal = diam_longitudinal_1**2 * math.pi / 4
area_total_separada = n_barras_i * area_longitudinal for n_barras_i in barras_camadas_1
Ast = sum(area_total_separada)


# area_longitudinal = diam_longitudinal_1**2 * math.pi / 4
# #area_longitudinal = [diametro**2 * math.pi / 4 for diametro in diam_longitudinal_2]
# area_total_separada = [n_barras_i * area_i for n_barras_i, area_i in zip(barras_camadas_2, area_longitudinal)]
# Ast = sum(area_total_separada)
    # area_longitudinal = 3.1415/4*self.diam_longitudinal**2
    # n_barras = sum(self.barras_camadas)
    # Ast = n_barras * area_longitudinal
    
n_camadas = len(barras_camadas_2)
      
dlinha = cobrimento + diam_transversal_2 + diam_longitudinal_2/2
    # dlinha = self.cobrimento + self.diam_transversal + self.diam_longitudinal/2
slinha = espacamento_s_2 - diam_transversal_2
    
hx = b_beam - 2*cobrimento - diam_transversal_2   # dc
hy = h_beam - 2*cobrimento - diam_transversal_2   # bc
    
wlinha_x = ((hx - diam_transversal_2 - 2*diam_longitudinal_2)-(npernas_yy-2)*diam_longitudinal_2) / (npernas_yy-1)
wlinha_y = ((hy - diam_transversal_2 - 2*diam_longitudinal_2)-(npernas_xx-2)*diam_longitudinal_2) / (npernas_xx-1)

num_wlinha_x = 2 * (npernas_yy - 1)  
num_wlinha_y = 2 * (npernas_xx - 1) 
    
rho_cc = Ast/(hx * hy)  # ratio of area of longitudinal reinforcement to area of core of section 
Ae = (hx * hy - num_wlinha_x*(wlinha_x**2)/6 - num_wlinha_y*(wlinha_y**2)/6)*(1-0.5*slinha/hx)*(1-0.5*slinha/hy)   # Area of effectively confined concrete core
Acc = (hx * hy) * (1 - rho_cc) # Area of core concrete without longitudinal reinforcement (le resta área al núcle de concreto confinado hx*hy ó dc*bc)
ke = Ae / Acc    # Coefinement effectiveness coefficient
    
rho_x = (npernas_xx * math.pi * diam_transversal_2**2 / 4) / (espacamento_s_2*hy)   # Total area of transverse reinforcement parallel to the x axis
rho_y = (npernas_yy * math.pi * diam_transversal_2**2 / 4) / (espacamento_s_2*hx)   # Total area of transverse reinforcement parallel to the y axis
    
flx = ke * rho_x * fy_transv_1 * (10**6) # Lateral confining stress on the concrete in the direction x
fly = ke * rho_y * fy_transv_1 * (10**6) # Lateral confining stress on the concrete in the direction y
q = np.minimum(flx,fly)/np.maximum(flx,fly) 
    
A = 6.8886 - (0.6069 + 17.275*q)*np.exp(-4.989*q)
B = ((4.5)/(5/A*(0.9849-0.6306*np.exp(-3.8939*q))-0.1))-5
xlinha = (flx + fly) / (2*(fc_inicial))
k1 = A*(0.1+0.9/(1+B*xlinha))
k2 = 5*k1  # para resistencia normal da armadura transversal (5k1)
    # k2 = 3*k1  # para elevada resistencia da armadura transversal (3k1)
fator_fc_confinado = 1+k1*xlinha
fator_ec_confinado = 1+k2*xlinha
fc_confinado_2 = fc_inicial*fator_fc_confinado          #Resistência máxima do concreto confinado à compressão
ec_confinado_2 = epsc*fator_ec_confinado         #Deformação de pico do concreto confinado à compressão 
if fc_confinado > fc_inicial:
    xn = 30                                                      #ok
    n = Ec*ec_confinado/(fc_confinado)                           #ok
    r = n/(n-1)                                                  #ok
else:
    fc_confinado = fc_inicial
    ec_confinado = epsc
    xn = 2.3
    


epsc_u_NC2 = 0.008  #0.008 La defomración unitaria de ruptura del hormigón CONFINADO no es dado en el paper 34. Se asume 0.008. 

    # ops.uniaxialMaterial('Concrete01', matTag, fpc, epsc0, fpcu, epsU)
ops.uniaxialMaterial('Concrete01', ConcretoC_Tag_2, -fc_confinado_2, -ec_confinado_2, -fc_confinado_2*0.20, -epsc_u_NC2)  


#CONFINADO 2 - 88cm
matTag_C = ConcretoC_Tag_2
fpc_c = -fc_confinado_2
epsc0_c = -ec_confinado_2
epsU_c = -epsc_u_NC2
ops.uniaxialMaterial('Concrete01', matTag_C, fpc_c, epsc0_c, fpcu_c, epsU_c)

    # #CONFINADO 2 - 88cm 
    # matTag_C = 4
    # fpc_c = 37.33e6 #N/m2
    # epsc0_c = 0.0026 #E0 = 25.6 GPa
    # fpcu_c = 17.85e6 #N/m2 aproximado do grafico 
    # epsU_c = 0.03 #do grafico
    # ops.uniaxialMaterial('Concrete01', matTag_C, fpc_c, epsc0_c, fpcu_c, epsU_c)

#Nós
ops.node(1, 0, 0)
ops.node(2, 0, 0.55)
ops.node(3, 0, 1.43)
ops.node(4, 0, 1.65)

ops.fix(1, 1, 1, 1)

#SEÇÃO FIBER - 5 fibras
c = 0.04 # cover (m)
# area_s = 0.0006446 
d = 0.02 #0.02 m de diam
area_s = np.pi*(d**2)/4 #barra de aço - longitudinal (m)

secTag = 1
fib_sec_1 = [['section', 'Fiber', 1, '-GJ', 1.0e6],
             ['patch', 'quad', 2, 4, 4, (-H/2)+c, (-B/2)+c, (H/2)-c, (-B/2)+c, (H/2)-c, (B/2)-c, (-H/2)+c, (B/2)-c], #confinado - 5
             ['patch', 'quad', 3, 2, 6, (-H/2), (-B/2), (-H/2)+c, (-B/2), (-H/2)+c, (B/2), (-H/2), (B/2)], #1 não confinado
             ['patch', 'quad', 3, 6, 2, (-H/2)+c, (B/2)-c, (H/2)-c, (B/2)-c, (H/2)-c, (B/2), (-H/2)+c, (B/2)], #2
             ['patch', 'quad', 3, 6, 2, (-H/2)+c, (-B/2), (H/2)-c, (-B/2), (H/2)-c, (-B/2)+c, (-H/2)+c, (-B/2)+c], #3
             ['patch', 'quad', 3, 2, 6, (H/2)-c, (-B/2), (H/2), (-B/2), (H/2), (B/2), (H/2)-c, (B/2)], #4
             ['layer', 'straight', 1, 4, area_s, (H/2)-c, (B/2)-c, (H/2)-c, (-B/2)+c],
             ['layer', 'straight', 1, 2, area_s, (H-2*c)/6, (B/2)-c, (H-2*c)/6, (-B/2)+c],
             ['layer', 'straight', 1, 2, area_s, -(H-2*c)/6, (B/2)-c, -(H-2*c)/6, (-B/2)+c],
             ['layer', 'straight', 1, 4, area_s, (-H/2)+c, (B/2)-c, (-H/2)+c, (-B/2)+c]]
opsv.fib_sec_list_to_cmds(fib_sec_1)
matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
opsv.plot_fiber_section(fib_sec_1, matcolor=matcolor)
plt.axis('equal')
plt.savefig('fibsec_rc.png')
plt.show()



#TRANSFORMAÇÃO
transfTag = 1
ops.geomTransf('Linear', transfTag)
# ops.geomTransf('PDelta', transfTag) considera a não linearidade geométrica do material

#INTEGRAÇÃO
integrationTag = 1       
ops.beamIntegration('Lobatto', integrationTag, secTag, 6)

#ELEMENTOS
eleTagAB = 1
eleTagBC = 2
eleTagCD = 3
ops.element('forceBeamColumn', eleTagAB, 1, 2, transfTag, integrationTag) #AB
ops.element('forceBeamColumn', eleTagBC, 2, 3, transfTag, integrationTag) #BC
ops.element('forceBeamColumn', eleTagCD, 3, 4, transfTag, integrationTag) #CD

#PESO PRÓPRIO
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)

rho = 2500
g = 9.81
A = B * H
peso_proprio = -rho * A * g * 0,4125 #L/n de nós

for i in range(2, 4 + 2):
    ops.load(i, 0.0, peso_proprio, 0.0)

#FORÇA AXIAL
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)

P = -(0.1*fpc*A) # kN vertical para baixo #conferir fpc 

ops.load(4, 0.0, P, 0.0)   # (Fx, Fy, Mz)

#TIMESERIES 
timeTag = 1
dt = 0.01
ops.timeSeries('Path', timeTag, '-filePath', 'forca.txt', '-dt', dt)
ops.pattern('Plain', 1, 1)
ops.load(4, 1.0, 0.0, 0.0)  # carga horizontal no topo

# ANÁLISE
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

with open('forca.txt', 'r') as f:
    lines = [float(line.strip()) for line in f]
nSteps = len(lines)

#ARMAZENAR
data = np.zeros((nSteps+1, 2))
data[0,:] = [0,0]

for j in range(nSteps):
    ok = ops.analyze(1, dt)
    if ok != 0:
        print(f"⚠️ Falha na análise no passo {j}")
        break

    disp = ops.nodeDisp(4, 1)
    force = lines[j]
    data[j+1,:] = [disp, force]

#SALVAR
df = pd.DataFrame(data, columns=['Deslocamento (m)', 'Força (N)'])
df.to_excel('resultados_pilar_feio.xlsx', index=False, engine='openpyxl')
print("✅ Arquivo 'resultados_pilar_feio.xlsx' criado com sucesso!")

#PLOTAR GRÁFICO
plt.plot(data[:,0], data[:,1], 'r-')
plt.xlabel('Deslocamento (m)')
plt.ylabel('Força (N)')
plt.title('Curva Força x Deslocamento')
plt.grid(True)
plt.show()

opsv.plot_model(element_labels=0, node_labels=0)
plt.show()
