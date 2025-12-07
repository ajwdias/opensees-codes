import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import pandas as pd 


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

#NAO CONFINADO curva 5 e 8
#uniaxialMaterial('Concrete01', matTag, fpc, epsc0, fpcu, epsU)
matTag_NC = 2
fpc = 32*(10**6) #Pa
epsc0 = 0.0025 #E0 = 25.6 GPa
fpcu = 0.2 * fpc #aproximado 
epsU = 0.0045 #aproximado do grafico
ops.uniaxialMaterial('Concrete01', matTag_NC, fpc, epsc0, fpcu, epsU)

#CONFINADO curva 5 e 6
matTag_C = 3
fpc_c = 37.33e6 #Mpa
epsc0_c = 0.0026 #E0 = 25.6 GPa
fpcu_c = 17.85e6 #aproximado N/m²
epsU_c = 0.08 #aproximado ??
ops.uniaxialMaterial('Concrete01', matTag_C, fpc_c, epsc0_c, fpcu_c, epsU_c)

#Nós
ops.node(1, 0, 0)
ops.node(2, 0, 0.55)
ops.node(3, 0, 1.10)
ops.node(4, 0, 1.65)

ops.fix(1, 1, 1, 1)

#SEÇÃO FIBER - 5 fibras
c = 0.05 # cover (m)
# area_s = 0.0006446 
d = 0.02 #20mm de diam
area_s = np.pi*(d**2)  #barra de aço - longitudinal (m)

secTag = 1
fib_sec_1 = [['section', 'Fiber', 1, '-GJ', 1.0e6],
             ['patch', 'quad', 2, 4, 4, (-H/2)+c, (-B/2)+c, (H/2)-c, (-B/2)+c, (H/2)-c, (B/2)-c, (-H/2)+c, (B/2)-c], #confinado - 5
             ['patch', 'quad', 3, 2, 6, (-H/2), (-B/2), (-H/2)+c, (-B/2), (-H/2)+c, (B/2), (-H/2), (B/2)], #1 não confinado
             ['patch', 'quad', 3, 6, 2, (-H/2)+c, (B/2)-c, (H/2)-c, (B/2)-c, (H/2)-c, (B/2), (-H/2)+c, (B/2)], #2
             ['patch', 'quad', 3, 6, 2, (-H/2)+c, (-B/2), (H/2)-c, (-B/2), (H/2)-c, (-B/2)+c, (-H/2)+c, (-B/2)+c], #3
             ['patch', 'quad', 3, 2, 6, (H/2)-c, (-B/2), (H/2), (-B/2), (H/2), (B/2), (H/2)-c, (B/2)], #4
             ['layer', 'straight', 1, 3, area_s, (H/2)-c, (B/2)-c, (H/2)-c, (-B/2)+c],
             ['layer', 'straight', 1, 2, area_s, 0, (B/2)-c, 0, (-B/2)+c],
             ['layer', 'straight', 1, 3, area_s, (-H/2)+c, (B/2)-c, (-H/2)+c, (-B/2)+c]]
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

#TIMESERIES 
timeTag = 1
dt = 0.01
ops.timeSeries('Path', timeTag, '-filePath', 'forca_tanaka.txt', '-dt', dt)
ops.pattern('Plain', 1, 1)
ops.load(4, 1.0, 0.0, 0.0)  # carga horizontal no topo

# ANÁLISE
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

with open('forca_tanaka.txt', 'r') as f:
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

    disp = ops.nodeDisp(2, 1)
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
