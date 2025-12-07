import openseespy.opensees as ops
import opsvis as vis
import matplotlib.pyplot as plt
import numpy as np

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)

# geometria
altura = 3.0
base = 0.3
n_ele = 10
delta_z = altura / n_ele

# materiais
ops.uniaxialMaterial("Concrete01", 1, -30e6, -0.002, -20e6, -0.006)
ops.uniaxialMaterial("Steel01", 2, 500e6, 200e9, 0.01)

# fiber
A_barra = 0.0003
ops.section("Fiber", 1)
ops.patch("rect", 1, 10, 10, -base/2, -base/2, base/2, base/2)
ops.layer("straight", 2, 4, A_barra, -0.1, -0.1, -0.1, 0.1)
ops.layer("straight", 2, 4, A_barra, 0.1, -0.1, 0.1, 0.1)

ops.beamIntegration("Lobatto", 1, 1, 5)
ops.geomTransf("Linear", 1)

# nós
for i in range(n_ele + 1):
    ops.node(i+1, 0.0, i*delta_z)
ops.fix(1, 1, 1, 1)

# peso próprio
peso_linear = base * base * 25e3
ops.timeSeries("Constant", 1)
ops.pattern("Plain", 1, 1)
for i in range(2, n_ele+2):
    ops.load(i, 0.0, -peso_linear * delta_z, 0.0)

# elementos
for i in range(n_ele):
    ops.element("forceBeamColumn", i+1, i+1, i+2, 1, 1)

# série temporal senoidal (c/path)
periodo = 2.0
dt = 0.05
n_steps = 200
tempo_total = dt * n_steps
amplitude = 1.0

tempos = np.arange(0, tempo_total + dt, dt)
valores = amplitude * np.sin(2 * np.pi * tempos / periodo)

ops.timeSeries("Path", 2, "-time", *tempos, "-values", *valores)

ops.pattern("Plain", 2, 2)
ops.load(n_ele+1, 10000.0, 0.0, 0.0)

# sistema e algoritmos
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.test("NormDispIncr", 1e-8, 10)
ops.algorithm("Newton")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")

# armazenar resultados
times = []
ux_top = []
uy_top = []

# loop de análise
for i in range(n_steps):
    ops.analyze(1, dt)
    disp = ops.nodeDisp(n_ele+1)
    times.append((i+1)*dt)
    ux_top.append(disp[0])
    uy_top.append(disp[1])

# deslocamento final
ux_final, uy_final, _ = ops.nodeDisp(n_ele+1)
print(f"Deslocamento no topo após {tempo_total:.2f} s:")
print(f"  Horizontal (ux): {ux_final:.6f} m")
print(f"  Vertical   (uy): {uy_final:.6f} m")

# gráfico horizontal
plt.figure()
plt.plot(times, ux_top)
plt.xlabel("Tempo [s]")
plt.ylabel("Deslocamento horizontal [m]")
plt.title("Deslocamento horizontal do topo")
plt.grid(True)
plt.show()

# gráfico vertical
plt.figure()
plt.plot(times, uy_top)
plt.xlabel("Tempo [s]")
plt.ylabel("Deslocamento vertical [m]")
plt.title("Deslocamento vertical do topo")
plt.grid(True)
plt.show()

# visualização da deformada final
vis.plot_defo()
