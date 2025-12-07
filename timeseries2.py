import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# 1️⃣ Modelo
# -----------------------------
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Parâmetros do material e geometria
matTag_aco = 1
B = 0.3
H = 0.3
I = (B*H**3)/12
L = 3
E = 210e9
k0 = (3*E*I)/L**3

# Material Steel01 (elasticidade + plastificação)
Fy = 5e5       # força de escoamento
b  = 0.01      # razão de pós-escoamento
ops.uniaxialMaterial('Steel01', matTag_aco, Fy, k0, b)

# Nós e elemento zeroLength
ops.node(1, 0, 0)
ops.node(2, 0, 0)
ops.fix(1, 1, 1, 1)
ops.fix(2, 0, 1, 1)
ops.element('zeroLength', 1, 1, 2, '-mat', matTag_aco, '-dir', 1)

# -----------------------------
# 2️⃣ TimeSeries cíclica
# -----------------------------
nSteps = 3000
nCycles = 3
Pmax = 1e6  # força máxima > Fy para gerar plastificação

t = np.linspace(0, nCycles, nSteps+1)
# força senoidal para gerar loops
val = Pmax * np.sin(2 * np.pi * nCycles * t)

ops.timeSeries('Path', 1, '-time', *t, '-values', *(val))
ops.pattern('Plain', 1, 1)
ops.load(2, 1.0, 0, 0)

# -----------------------------
# 3️⃣ Configuração da análise
# -----------------------------
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.test('NormDispIncr', 1e-6, 200)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

dt = t[1] - t[0]

# -----------------------------
# 4️⃣ Execução da análise
# -----------------------------
data = np.zeros((nSteps+1, 2))
data[0, :] = [0.0, 0.0]

for j in range(nSteps):
    ok = ops.analyze(1, dt)
    if ok != 0:
        print(f"⚠️ Falha no passo {j}")
        break
    disp = ops.nodeDisp(2, 1)
    force = ops.eleForce(1)[0]
    data[j+1, :] = [disp, force]

# -----------------------------
# 5️⃣ Plot da curva histerética
# -----------------------------
plt.figure(figsize=(6,5))
plt.plot(data[:,0], data[:,1], '-b')
plt.xlabel('Deslocamento (m)')
plt.ylabel('Força (N)')
plt.title(f'Curva Histerética - {nCycles} ciclos')
plt.grid(True)
plt.show()
