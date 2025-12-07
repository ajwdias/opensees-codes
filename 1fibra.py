import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import opsvis as opsv

# Modelo
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

#PARÂMETROS
matTag_aco = 1
B = 0.3 #m
H = 0.3 #m
A = B*H
I = (B*(H**3))/12 #m4
L_col = 3 #m
H_viga = 3 #m
Fy = 500e6   # N/m²
E = 210e9 #Pa
b = 0.003226
ops.uniaxialMaterial('Steel01', matTag_aco, Fy, E, b)

#NÓS
ops.node(1, 0, 0) #A
ops.node(2, 0, H_viga) #B
ops.node(3, L_col, H_viga) #C
ops.node(4, L_col, 0) #D

ops.fix(1, 1, 1, 0)
ops.fix(4, 1, 1, 0)

#SEÇÃO FIBER - 1 fibra
secTag = 1
fib_sec_1 = [['section', 'Fiber', 1, '-GJ', 1.0e6],
    ['patch', 'quad', 1, 4, 4, -H/2, -B/2, H/2, -B/2, H/2, B/2, -H/2, B/2], 
             ]
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
ops.beamIntegration('Lobatto', integrationTag, secTag, 5)

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
ops.timeSeries('Path', timeTag, '-filePath', 'forca.txt', '-dt', dt)
ops.pattern('Plain', 1, 1)
ops.load(2, 1.0, 0.0, 0.0)  # carga horizontal no topo

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

    disp = ops.nodeDisp(2, 1)
    force = lines[j]
    data[j+1,:] = [disp, force]

#SALVAR
df = pd.DataFrame(data, columns=['Deslocamento (m)', 'Força (N)'])
df.to_excel('resultados_portico.xlsx', index=False, engine='openpyxl')
print("✅ Arquivo 'resultados_portico.xlsx' criado com sucesso!")

#PLOTAR GRÁFICO
plt.plot(data[:,0], data[:,1], 'r-')
plt.xlabel('Deslocamento (m)')
plt.ylabel('Força (N)')
plt.title('Curva Força x Deslocamento')
plt.grid(True)
plt.show()

opsv.plot_model(element_labels=0, node_labels=0)
plt.show()
