import openseespy.opensees as ops
import opsvis as vis
import matplotlib.pyplot as plt
import math

ops.wipe()

ops.model('basic', '-ndm', 2, '-ndf', 3)

# Materiais
# Concrete01: f'c = -30 MPa, eps_c0 = -0.002, f'cu = -20 MPa, eps_u = -0.006
ops.uniaxialMaterial("Concrete01", 1, -30e6, -0.002, -20e6, -0.006)

# Steel01: Fy=500 MPa, E=200 GPa, b=0.01 (hardening ratio)
ops.uniaxialMaterial("Steel01", 2, 500e6, 200e9, 0.01)

# Geometria e discretização
b = 0.3    # largura seção (m)
h = 0.3    # altura seção (m)
cover = 0.04  # cobrimento (m)
n_fibers_y = 10
n_fibers_z = 10

# Criar seção fiber
sectionTag = 1
ops.section("Fiber", sectionTag)

# Concrete patch (região interna, descontado cobrimento)
ops.patch("rect", 1, n_fibers_y, n_fibers_z, 
          -b/2 + cover, -h/2 + cover, b/2 - cover, h/2 - cover)

# Barras de aço
A_bar = math.pi * (0.016**2) / 4  # área barra de 16mm de diâmetro

# Coordenadas das barras (4 cantos)
bar_coords = [
    (-b/2 + cover, -h/2 + cover),
    (b/2 - cover, -h/2 + cover),
    (-b/2 + cover, h/2 - cover),
    (b/2 - cover, h/2 - cover)
]

# Camada de barras com aço Steel01 (material tag 2)
ops.layer("straight", 2, len(bar_coords), A_bar, 
          *[coord for bar in bar_coords for coord in bar])

# Criar nós
H = 3.0
n_ele = 10
delta_z = H / n_ele

for i in range(n_ele + 1):
    ops.node(i+1, 0.0, i*delta_z)

# Fixar base
ops.fix(1, 1, 1, 1)

# Transformação geométrica
ops.geomTransf("Linear", 1)

# Definir regra de integração
ops.beamIntegration("Lobatto", 1, sectionTag, 5)

# Criar elementos forceBeamColumn com regra de integração 1
for i in range(n_ele):
    ops.element("forceBeamColumn", i+1, i+1, i+2, sectionTag, 1)


# Cargas
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)

rho = 2500
g = 9.81
A = b * h
peso_proprio = -rho * A * g * delta_z

for i in range(2, n_ele + 2):
    ops.load(i, 0.0, peso_proprio, 0.0)

# Força horizontal no topo
ops.load(n_ele + 1, 10000.0, 0.0, 0.0)

# Análise
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Newton")
ops.analysis("Static")

ops.analyze(1)

# Resultados
ux = ops.nodeDisp(n_ele+1, 1)
uy = ops.nodeDisp(n_ele+1, 2)
print(f"Deslocamento topo: ux = {ux:.6f} m, uy = {uy:.6f} m")

# Visualização 
vis.plot_model(node_labels=False, element_labels=False)
vis.plot_defo()

# Gráfico deformada matplotlibS
x_def = []
y_def = []

for i in range(n_ele + 1):
    x = ops.nodeCoord(i+1)[0] + ops.nodeDisp(i+1)[0]
    y = ops.nodeCoord(i+1)[1] + ops.nodeDisp(i+1)[1]
    x_def.append(x)
    y_def.append(y)

plt.figure()
plt.plot(x_def, y_def, 'b-o')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('Deformada do Pilar (Fibra)')
plt.grid(True)
plt.axis('equal')
plt.show()


