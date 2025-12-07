import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

# ------------------------
# Limpar e criar modelo
# ------------------------
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# ------------------------
# PARÂMETROS (nomes claros)
# ------------------------
matTag_aco = 1
B = 0.3           # largura (m)
H = 0.3           # altura (m)
A = B * H
Izz = (B * (H**3)) / 12.0
L_col = 3.0       # vão em x (m)
H_viga = 3.0      # altura das colunas (m)
Fy = 500e6        # Pa (tensão de escoamento)
E = 210e9         # Pa
b_hard = 0.01     # coeficiente encruamento (Steel01)

# Material
ops.uniaxialMaterial('Steel01', matTag_aco, Fy, E, b_hard)

# ------------------------
# Nós
# ------------------------
ops.node(1, 0.0, 0.0)           # A (base esquerda)
ops.node(2, 0.0, H_viga)        # B (topo esquerda)
ops.node(3, L_col, H_viga)      # C (topo direita)
ops.node(4, L_col, 0.0)         # D (base direita)

# Apoios: engastes nas bases (transl + rot travados)
ops.fix(1, 1, 1, 1)
ops.fix(4, 1, 1, 1)

# ------------------------
# SEÇÃO FIBER
# ------------------------
secTag = 1
ops.section('Fiber', secTag)
nY = 10
nZ = 10
y_coords = np.linspace(-H/2, H/2, nY)
z_coords = np.linspace(-B/2, B/2, nZ)
A_fibra = (H / nY) * (B / nZ)
for y in y_coords:
    for z in z_coords:
        ops.fiber(y, z, A_fibra, matTag_aco)

# ------------------------
# Transform e integração
# ------------------------
transfTag = 1
ops.geomTransf('Linear', transfTag)

integrationTag = 1
nIntPts = 5
ops.beamIntegration('Lobatto', integrationTag, secTag, nIntPts)

# ------------------------
# Elementos forceBeamColumn
# ------------------------
ops.element('forceBeamColumn', 1, 1, 2, transfTag, integrationTag)  # AB
ops.element('forceBeamColumn', 2, 2, 3, transfTag, integrationTag)  # BC
ops.element('forceBeamColumn', 3, 3, 4, transfTag, integrationTag)  # CD

# ------------------------
# Ler arquivo de forças (robusto)
# ------------------------
force_file = 'forca.txt'
loads = []
with open(force_file, 'r') as f:
    for line in f:
        s = line.strip()
        if s == '' or s.startswith('#'):
            continue
        try:
            val = float(s)
            loads.append(val)
        except ValueError:
            print(f"⚠️ Linha ignorada no arquivo de forças (não numérica): '{s}'")

if len(loads) == 0:
    raise RuntimeError("Arquivo de forças vazio ou sem linhas válidas. Verifique 'forca.txt'.")

loads = np.array(loads, dtype=float)
nSteps = len(loads)

# ------------------------
# TimeSeries / Pattern
# ------------------------
timeTag = 1
dt = 0.01  # tempo entre valores no arquivo (s)
ops.timeSeries('Path', timeTag, '-dt', dt, '-values', *loads)
ops.pattern('Plain', 1, timeTag)
# aplica na posição do nó B (topo esquerdo) ou escolha outro nó
ops.load(2, 1.0, 0.0, 0.0)

# ------------------------
# Configuração da análise
# ------------------------
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')

# Adicionar pequeno amortecimento Rayleigh (ajustável)
# ops.rayleigh(alphaM, betaKcurr, betaKinit, betaKcomm)
ops.rayleigh(0.0, 1e-6, 0.0, 0.0)  # betaKcurr pequeno ajuda estabilidade

# Use Newmark (implicit) para integração direta
beta = 0.25
gamma = 0.5

ops.algorithm('Newton')    # algoritmo principal
ops.integrator('Newmark', gamma, beta)
ops.analysis('Transient')

# ------------------------
# Função utilitária para tentar analisar com fallback
# ------------------------
def try_analyze(tstep, dt_step, n_substeps=1):
    """
    Tenta executar ops.analyze(n_substeps, dt_step) e faz fallback:
    - primeiro tenta com algoritmo atual
    - se falhar, tenta ModifiedNewton
    - se ainda falhar, subdivide o passo em n_substeps (maior) subpassos
    Retorna True se ok, False se falhou definitivamente.
    """
    # tentativa direta
    ok = ops.analyze(tstep, dt_step)
    if ok == 0:
        return True

    # tentar ModifiedNewton once
    ops.algorithm('ModifiedNewton')
    ok = ops.analyze(tstep, dt_step)
    if ok == 0:
        # restaurar algoritmo para Newton para próximos passos
        ops.algorithm('Newton')
        return True

    # tentar subdividir em passos menores
    ops.algorithm('ModifiedNewton')
    sub = max(2, n_substeps)
    small_dt = dt_step / sub
    for k in range(sub):
        ok = ops.analyze(1, small_dt)
        if ok != 0:
            # falhou em subdivisão
            return False
    # se chegar aqui, ok
    ops.algorithm('Newton')
    return True

# ------------------------
# Loop de análise (com salvamento)
# ------------------------
data = np.zeros((nSteps+1, 2))
data[0, :] = [ops.nodeDisp(2, 1), 0.0]  # deslocamento inicial no nó 2 e força 0

failed_at = None
for j in range(nSteps):
    # tentamos avançar 1 passo de tempo dt
    success = try_analyze(1, dt, n_substeps=8)
    if not success:
        print(f"⚠️ Falha na análise no passo {j} (tempo ~ {(j+1)*dt:.2f} s). Tentativas de recuperação esgotadas.")
        failed_at = j
        break

    # colher resultados (deslocamento horizontal no nó 2)
    disp = ops.nodeDisp(2, 1)
    force = loads[j]
    data[j+1, :] = [disp, force]

# ------------------------
# Salvar resultados (xlsx se openpyxl disponível, senão csv)
# ------------------------
df = pd.DataFrame(data[: (0 if failed_at is None else failed_at+2) if False else (failed_at+2 if failed_at is not None else nSteps+1)],
                  columns=['Deslocamento (m)', 'Força (N)'])

# tentativa de salvar em excel, fallback para csv
try:
    df.to_excel('resultados_portico.xlsx', index=False, engine='openpyxl')
    print("✅ Arquivo 'resultados_portico.xlsx' criado com sucesso!")
except Exception as e:
    print(f"⚠️ Não foi possível salvar .xlsx ({e}). Salvando CSV em 'resultados_portico.csv'...")
    df.to_csv('resultados_portico.csv', index=False)

# ------------------------
# Plot
# ------------------------
plt.figure(figsize=(7,5))
plt.plot(data[:df.shape[0], 0], data[:df.shape[0], 1], '-o', markersize=3)
plt.xlabel('Deslocamento (m)')
plt.ylabel('Força (N)')
plt.title('Curva Força x Deslocamento - Coluna de Aço Não Linear')
plt.grid(True)
plt.show()
