import openseespy.opensees as ops
import opsvis as vis
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------
# MODELO
# ----------------------------
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Parâmetros do pilar
h_pilar = 3.0  # altura do pilar (m)
nós = 10
dx = h_pilar / nós
b = 0.3        # largura da seção (m)
h = 0.5        # altura da seção (m)
cover = 0.03   # cobrimento (m)

# ----------------------------
# MATERIAIS UNIAXIAIS
# ----------------------------
ops.uniaxialMaterial('Concrete02', 1, -41.3e6, -0.002, 0.2*-41.3e6, -0.006)
ops.uniaxialMaterial('Concrete02', 2, -43558991.97, 0.00254697, -8260000.0, -0.008)
ops.uniaxialMaterial('Steel01', 3, 500e6, 200e9, 0.01)

# ----------------------------
# SEÇÃO FIBER
# ----------------------------
ops.section('Fiber', 1)
# Cobertura de concreto
ops.patch('rect', 1, 10, 1, -b/2, h/2, b/2, h/2 - cover)   # topo
ops.patch('rect', 1, 10, 1, -b/2, -h/2 + cover, b/2, -h/2) # base
ops.patch('rect', 1, 1, 10, -b/2, -h/2 + cover, -b/2 + cover, h/2 - cover)  # esquerda
ops.patch('rect', 1, 1, 10, b/2 - cover, -h/2 + cover, b/2, h/2 - cover)    # direita
# Núcleo confinado
ops.patch('rect', 2, 10, 10, -b/2 + cover, -h/2 + cover, b/2 - cover, h/2 - cover)
# Barras de aço
num_barras = 2
area_barra = np.pi * (0.012)**2 / 4
ops.layer('straight', 3, num_barras, area_barra, -b/2 + cover, h/2 - cover, b/2 - cover, h/2 - cover) # topo
ops.layer('straight', 3, num_barras, area_barra, -b/2 + cover, -h/2 + cover, b/2 - cover, -h/2 + cover) # base

# ----------------------------
# TRANSFORMAÇÃO E INTEGRAÇÃO
# ----------------------------
ops.geomTransf("Linear", 1)
ops.beamIntegration("Lobatto", 1, 1, 3)

# ----------------------------
# NÓS E APOIO
# ----------------------------
for i in range(nós + 1):
    ops.node(i + 1, 0.0, i * dx)
ops.fix(1, 1, 1, 1)

# ----------------------------
# ELEMENTOS
# ----------------------------
for i in range(nós):
    ops.element("forceBeamColumn", i + 1, i + 1, i + 2, 1, 1)

# ----------------------------
# MASSA DISTRIBUÍDA (necessário para análise dinâmica)
# ----------------------------
massa_linear = b * h * 25e3  # kg/m
for i in range(2, nós + 2):
    ops.mass(i, 0.0, massa_linear * dx, 0.0)

# ----------------------------
# PESO PRÓPRIO
# ----------------------------
peso_linear = b * h * 25e3  # N/m
ops.timeSeries("Constant", 1)
ops.pattern("Plain", 1, 1)
for i in range(2, nós + 2):
    ops.load(i, 0.0, -peso_linear * dx, 0.0)

# ----------------------------
# CARGA DINÂMICA SENOIDAL NO TOPO
# ----------------------------
periodo = 2.0
dt = 0.001       # passo de tempo menor para regime plástico
n_steps = 5000
amplitude = 21000  # N

tempos = np.arange(0, dt*n_steps + dt, dt)
valores = amplitude * np.sin(2 * np.pi * tempos / periodo)

ops.timeSeries("Path", 2, "-time", *tempos, "-values", *valores)
ops.pattern("Plain", 2, 2)
ops.load(nós + 1, 1.0, 0.0, 0.0)  # força multiplicada pelo timeSeries

# ----------------------------
# CONFIGURAÇÃO DA ANÁLISE
# ----------------------------
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
# integrador HHT para estabilidade
ops.integrator("HHT", -0.05, 0.25, 0.5)

# ----------------------------
# LOOP DE ANÁLISE ROBUSTO
# ----------------------------
dt_current = dt
min_dt = 1e-5
reduction_factor = 0.5

# fallback em cadeia
def try_analyze_fallback(nSteps, dt_step):
    strategies = [
        ("BandGeneral", ("NormDispIncr", 1e-8, 10), "Newton"),
        ("BandGeneral", ("NormDispIncr", 1e-6, 20), "ModifiedNewton"),
        ("BandGeneral", ("NormDispIncr", 1e-6, 50), "NewtonLineSearch"),
        ("BandGeneral", ("NormDispIncr", 1e-5, 100), "BFGS"),
    ]
    for systemType, testArgs, algoType in strategies:
        ops.system(systemType)
        ops.test(*testArgs)
        ops.algorithm(algoType)
        ok = ops.analyze(nSteps, dt_step)
        if ok == 0:
            return 0
    return -1

# Armazenamento de resultados
times = []
ux_top = []
uy_top = []
deslocamentos = []
forcas = []

i = 0
while i < n_steps:
    ok = try_analyze_fallback(1, dt_current)

    if ok == 0:
        # registro de resultados válidos
        disp = ops.nodeDisp(nós + 1)
        ux = disp[0]
        uy = disp[1]

        times.append(i * dt_current)
        ux_top.append(ux)
        uy_top.append(uy)
        deslocamentos.append(ux)

        ops.reactions()
        rx = ops.nodeReaction(1)[0]
        forcas.append(rx)

        i += 1  # avançar passo
    else:
        # reduzir dt se falhou
        if dt_current > min_dt:
            dt_current *= reduction_factor
            print(f"Falha no passo {i}, reduzindo dt para {dt_current:.6f} s e tentando novamente...")
        else:
            print(f"Falha crítica no passo {i}, dt já no mínimo {min_dt}. Interrompendo análise.")
            break

# ----------------------------
# RESULTADOS
# ----------------------------
ux_final, uy_final, _ = ops.nodeDisp(nós + 1)
print(f"Deslocamento final topo: ux = {ux_final:.6f} m, uy = {uy_final:.6f} m")

# ----------------------------
# GRÁFICOS
# ----------------------------
plt.figure()
plt.plot(times, ux_top)
plt.xlabel("Tempo [s]")
plt.ylabel("Deslocamento horizontal [m]")
plt.title("Deslocamento horizontal do topo")
plt.grid(True)

plt.figure()
plt.plot(times, uy_top)
plt.xlabel("Tempo [s]")
plt.ylabel("Deslocamento vertical [m]")
plt.title("Deslocamento vertical do topo")
plt.grid(True)

plt.figure()
plt.plot(deslocamentos, forcas, '-o', markersize=2)
plt.xlabel("Deslocamento horizontal no topo [m]")
plt.ylabel("Força horizontal na base [N]")
plt.title("Curva Força × Deslocamento")
plt.grid(True)

plt.show()

# ----------------------------
# VISUALIZAÇÃO DA DEFORMADA FINAL
# ----------------------------
vis.plot_defo()
