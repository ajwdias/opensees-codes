import openseespy.opensees as ops
import opsvis as vis
import matplotlib.pyplot as plt

# --- 1. Limpar modelo ---
ops.wipe()

# --- 2. Definir modelo (2D, 3 dof) ---
ops.model('basic', '-ndm', 2, '-ndf', 3)

# --- 3. Parâmetros geométricos e materiais ---
H = 3.0          # altura (m)
b = 0.3          # largura seção (m)
h = 0.3          # altura seção (m)
n_ele = 10       # número de elementos
E = 25e9         # módulo elasticidade (Pa)
rho = 2500       # densidade (kg/m³)
g = 9.81         # gravidade (m/s²)

A = b * h
Iz = (b * h**3) / 12
delta_z = H / n_ele

# --- 4. Criar nós e massa ---
for i in range(n_ele + 1):
    ops.node(i+1, 0.0, i * delta_z)
    ops.mass(i+1, 0.0, rho * A * delta_z, 0.0)  # massa no y (vertical)

# --- 5. Fixar base ---
ops.fix(1, 1, 1, 1)

# --- 6. Definir transformação geométrica ---
ops.geomTransf("Linear", 1)

# --- 7. Criar elementos elasticoBeamColumn ---
for i in range(n_ele):
    ops.element("elasticBeamColumn", i+1, i+1, i+2, A, E, Iz, 1)

# --- 8. Definir cargas ---
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)

peso_proprio = -rho * A * g * delta_z  # carga peso por elemento (negativo no y)

# Aplica peso próprio nos nós (exceto base)
for i in range(2, n_ele + 2):
    ops.load(i, 0.0, peso_proprio, 0.0)

# Força horizontal concentrada no topo
ops.load(n_ele + 1, 10000.0, 0.0, 0.0)

# --- 9. Configurar análise ---
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Newton")
ops.analysis("Static")

# --- 10. Executar análise ---
ops.analyze(1)

# --- 11. Resultados ---
ux = ops.nodeDisp(n_ele+1, 1)
uy = ops.nodeDisp(n_ele+1, 2)
print(f"Deslocamento no topo: ux = {ux:.6f} m, uy = {uy:.6f} m")

# --- 12. Visualização com opsvis ---
vis.plot_model(node_labels=False, element_labels=False)
vis.plot_defo()

# --- 13. Gráfico deformada com matplotlib ---
x_def = []
y_def = []

for i in range(n_ele + 1):
    x = ops.nodeCoord(i+1)[0] + ops.nodeDisp(i+1)[0]
    y = ops.nodeCoord(i+1)[1] + ops.nodeDisp(i+1)[1]
    x_def.append(x)
    y_def.append(y)

plt.figure()
plt.plot(x_def, y_def, 'b-o', label='Deformada')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('Deformada do Pilar (matplotlib)')
plt.grid(True)
plt.axis('equal')
plt.legend()
plt.show()
