import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

ops.wipe()

# ---------- modelo (2D, truss)
ops.model('basic', '-ndm', 2, '-ndf', 2)

# ---------- material (unidades físicas: Pa)
E = 210e9            # Pa
Fy = 250e6           # Pa (250 MPa)
b = 0.01             # razão pós-escoamento (adimensional)
# opcional: a1..a4 (isotrópico) se suportado
# use sem a1..a4 ou inclua se sua build aceitar
matTag = 1
ops.uniaxialMaterial('Steel01', matTag, Fy, E, b)  # Fy [Pa], E0 [Pa]

# ---------- geometria e seção (truss)
L = 1.0              # comprimento do elemento [m]
A = 1.0e-4           # área da seção [m^2] (100 mm^2)
eleTag = 1
ops.node(1, 0.0, 0.0)
ops.node(2, L, 0.0)
ops.fix(1, 1, 1)
ops.fix(2, 0, 1)     # nódulo livre em X, restrito Y

# truss: elemento axial que converte tensão(do material) em força = stress * A
ops.element('Truss', eleTag, 1, 2, A, matTag)

# ---------- informações úteis (para checar escala)
Fy_force = Fy * A                # força de escoamento [N]
k_axial = E * A / L              # rigidez axial [N/m]
u_yield = Fy_force / k_axial     # deslocamento de escoamento aproximado [m]
print(f"Fy (stress) = {Fy:.3e} Pa")
print(f"Área A = {A:.3e} m^2  -> força de escoamento Fy*A = {Fy_force:.3e} N")
print(f"Rigidez k = E*A/L = {k_axial:.3e} N/m")
print(f"Deslocamento de escoamento aproximado u_y = {u_yield:.6e} m ({u_yield*1e3:.3f} mm)")

# ---------- carregamento cíclico (força aplicada no nó 2, direção X)
nCycles = 4
nSteps = 1200
t = np.linspace(0, nCycles, nSteps)
# escolha Pmax > Fy*A para garantir plastificação
Pmax = 3.0 * Fy_force
print(f"Pmax aplicado = {Pmax:.3e} N  (Pmax/Fy*A = {Pmax/Fy_force:.2f})")

force = Pmax * np.sin(2*np.pi*t)    # vetor de forças [N]
np.savetxt('forca_truss.txt', force, fmt='%.6f')

# TimeSeries + Pattern (usa a série para escalar a carga definida abaixo)
ops.timeSeries('Path', 1, '-filePath', 'forca_truss.txt', '-dt', t[1]-t[0])
ops.pattern('Plain', 1, 1)
# define carga de referência 1.0 N no grau de liberdade X do nó 2 (será escalada pela timeSeries)
ops.load(2, 1.0, 0.0)

# ---------- análise estática incremental (quasi-estática)
ops.system('BandGeneral')
ops.numberer('Plain')
ops.constraints('Plain')
ops.test('NormUnbalance', 1e-6, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0)
ops.analysis('Static')

# ---------- loop com fallback
data = np.zeros((nSteps+1, 2))
data[0,:] = [0.0, 0.0]

for j in range(nSteps):
    ops.setTime(t[j])
    ok = ops.analyze(1)
    if ok != 0:
        # tentativas de recuperação (simples)
        print(f"⚠️ Falha no passo {j}, tentando NewtonLineSearch...")
        ops.algorithm('NewtonLineSearch', 0.8)
        ops.test('NormUnbalance', 1e-4, 100)
        ok = ops.analyze(1)

        if ok != 0:
            print(f"❌ Não convergiu no passo {j}. Parando análise.")
            break

    # deslocamento axial do nó 2 em X e força aplicada (valor do arquivo)
    disp = ops.nodeDisp(2, 1)    # deslocamento em X
    applied_force = force[j]
    data[j+1, :] = [disp, applied_force]

# ---------- salvar e plotar
df = pd.DataFrame(data, columns=['Deslocamento (m)', 'Força (N)'])
df.to_excel('resultados_truss.xlsx', index=False, engine='openpyxl')
print("Resultados salvos em 'resultados_truss.xlsx'.")

plt.figure(figsize=(7,5))
plt.plot(data[:,0], data[:,1], '-', lw=1)
plt.xlabel('Deslocamento (m)')
plt.ylabel('Força (N)')
plt.title('Histerese - Truss axial com Steel01 (unidades físicas)')
plt.grid(True)
plt.show()
