import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Modelo
# ----------------------------
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)  # 2D frame, 3 DOF/nó

# ----------------------------
# Nós
# ----------------------------
H = 0  # altura do pilar
ops.node(1, 0.0, 0.0)  # base
ops.node(2, 0.0, 0.0)    # topo
ops.fix(1, 1, 1, 1)    # base totalmente engastada

# ----------------------------
# Material Steel01
# ----------------------------
Fy = 36.0
E = 29000.0
b = 0.01
ops.uniaxialMaterial('Steel01', 1, Fy, E, b)

# ----------------------------
# Seção simplificada
# ----------------------------
A = 4.0      # área da seção
I = 1.0      # momento de inércia (suficiente para estabilidade)
ops.section('Fiber', 1)
ops.patch('rect', 1, 1, 1, -0.5, -0.5, 0.5, 0.5)

# ----------------------------
# Geometria e integração
# ----------------------------
ops.geomTransf('Linear', 1)
ops.beamIntegration('Lobatto', 1, 1, 5)  # 5 pontos de integração

# ----------------------------
# Elemento dispBeamColumn
# ----------------------------
ops.element('dispBeamColumn', 1, 1, 2, 1, 1)

# ----------------------------
# TimeSeries e LoadPattern
# ----------------------------
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 0.001, 0.0, 0.0)  # pequena referência para o integrador

# ----------------------------
# Sistema e análise
# ----------------------------
ops.system('ProfileSPD')
ops.numberer('Plain')
ops.constraints('Plain')
ops.algorithm('Newton')
ops.test('NormUnbalance', 1e-8, 10)
ops.analysis('Static')

# ----------------------------
# DisplacementControl no topo (horizontal)
# ----------------------------
Px = 160.0        # força final de referência
Nsteps = 200
ops.integrator('DisplacementControl', 2, 1, Px/Nsteps)  # nó 2, DOF horizontal

# ----------------------------
# Loop de análise
# ----------------------------
data = np.zeros((Nsteps+1, 2))
for i in range(Nsteps):
    ok = ops.analyze(1)
    if ok != 0:
        print(f'Convergência falhou no passo {i+1}')
        break
    ux = ops.nodeDisp(2, 1)
    Fx = ops.eleResponse(1, 'force')[0]
    data[i+1, 0] = ux
    data[i+1, 1] = Fx

# ----------------------------
# Gráfico
# ----------------------------
plt.figure(figsize=(6,4))
plt.plot(data[:, 0], data[:, 1], '-o', markersize=3)
plt.xlabel('Horizontal Displacement [mm]')
plt.ylabel('Horizontal Load [N]')
plt.title('Pilar de Aço Steel01 - Curva F–u')
plt.grid(True)
plt.show()
