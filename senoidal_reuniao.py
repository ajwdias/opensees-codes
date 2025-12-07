import openseespy.opensees as ops
import opsvis as vis
import matplotlib.pyplot as plt
import numpy as np

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Parâmetros do pilar
h_pilar = 3.0  # m
nós = 10
dx = h_pilar / nós
b = 0.3        # largura da seção (m)
h = 0.3        # altura da seção (m)
cover = 0.025  # cobrimento (m)

# MATERIAIS UNIAXIAIS

ops.uniaxialMaterial('Concrete02', 1, -36e6, -0.002, -20e6, -0.006)   # não confinado
ops.uniaxialMaterial('Concrete02', 2, -35e6, -0.0025, -25e6, -0.01)   # confinado
ops.uniaxialMaterial('Steel01', 3, 500e6, 200e9, 0.01)                # aço


# SEÇÃO FIBER

ops.section('Fiber', 1)

# Cobertura de concreto (não confinado)
ops.patch('rect', 1, 10, 1, -b/2, h/2, b/2, h/2 - cover)               # Topo
ops.patch('rect', 1, 10, 1, -b/2, -h/2 + cover, b/2, -h/2)             # Base
ops.patch('rect', 1, 1, 10, -b/2, -h/2 + cover, -b/2 + cover, h/2 - cover)  # Esquerda
ops.patch('rect', 1, 1, 10, b/2 - cover, -h/2 + cover, b/2, h/2 - cover)    # Direita

# Núcleo de concreto confinado
ops.patch('rect', 2, 10, 10, -b/2 + cover, -h/2 + cover, b/2 - cover, h/2 - cover)

# Barras de aço
num_barras = 2
area_barra = np.pi * (0.012)**2 / 4  # bitola 12 mm
ops.layer('straight', 3, num_barras, area_barra, -b/2 + cover, h/2 - cover,  b/2 - cover, h/2 - cover)  # topo
ops.layer('straight', 3, num_barras, area_barra, -b/2 + cover, -h/2 + cover, b/2 - cover, -h/2 + cover) # base

# TRANSFORMAÇÃO E INTEGRAÇÃO
ops.geomTransf("Linear", 1)
ops.beamIntegration("Lobatto", 1, 1, 3)

# NÓS E APOIO

for i in range(nós + 1):
    ops.node(i + 1, 0.0, i * dx)
ops.fix(1, 1, 1, 1)


# PESO PRÓPRIO
peso_linear = b * b * 25e3  # N/m
ops.timeSeries("Constant", 1)
ops.pattern("Plain", 1, 1)
for i in range(2, nós + 2):
    ops.load(i, 0.0, -peso_linear * dx, 0.0)


# ELEMENTOS
for i in range(nós):
    ops.element("forceBeamColumn", i + 1, i + 1, i + 2, 1, 1)


# CARGA SENOIDAL NO TOPO
periodo = 2.0
dt = 0.05
n_steps = 200
tempo_total = dt * n_steps
amplitude = 1.0

tempos = np.arange(0, tempo_total + dt, dt)
valores = amplitude * np.sin(2 * np.pi * tempos / periodo)

ops.timeSeries("Path", 2, "-time", *tempos, "-values", *valores)
ops.pattern("Plain", 2, 2)
ops.load(nós + 1, 2000.0, 0.0, 0.0)


# CONFIGURAÇÃO DA ANÁLISE
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.test("NormDispIncr", 1e-8, 10)
ops.algorithm("Newton")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")


# LOOP DE ANÁLISE (com força × deslocamento)
times = []
ux_top = []
uy_top = []
deslocamentos = []
forcas = []

for i in range(n_steps):
    ok = ops.analyze(1, dt)
    if ok != 0:
        print(f"Falha na análise no passo {i}")
        break
    
    disp = ops.nodeDisp(nós + 1)
    ux = disp[0]
    uy = disp[1]
    
    # Guardar deslocamentos por tempo
    times.append((i + 1) * dt)
    ux_top.append(ux)
    uy_top.append(uy)
    
    # Guardar deslocamento × força
    deslocamentos.append(ux)
    ops.reactions()
    rx = ops.nodeReaction(1)[0]  # reação horizontal na base
    forcas.append(rx)


# RESULTADOS NUMÉRICOS
ux_final, uy_final, _ = ops.nodeDisp(nós + 1)
print(f"Deslocamento no topo após {tempo_total:.2f} s:")
print(f"  Horizontal (ux): {ux_final:.6f} m")
print(f"  Vertical   (uy): {uy_final:.6f} m")


# GRÁFICOS
# Deslocamento horizontal × tempo
plt.figure()
plt.plot(times, ux_top)
plt.xlabel("Tempo [s]")
plt.ylabel("Deslocamento horizontal [m]")
plt.title("Deslocamento horizontal do topo")
plt.grid(True)

# Deslocamento vertical × tempo
plt.figure()
plt.plot(times, uy_top)
plt.xlabel("Tempo [s]")
plt.ylabel("Deslocamento vertical [m]")
plt.title("Deslocamento vertical do topo")
plt.grid(True)

# Força × deslocamento (histerese)
plt.figure()
plt.plot(deslocamentos, forcas, '-o', markersize=2)
plt.xlabel("Deslocamento horizontal no topo [m]")
plt.ylabel("Força horizontal na base [N]")
plt.title("Curva Força × Deslocamento")
plt.grid(True)

plt.show()


# VISUALIZAÇÃO DA DEFORMADA FINAL
vis.plot_defo()
