import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import math as math

# --- INICIALIZAÇÃO DO MODELO ---
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# --- 0. FUNÇÃO DE CÁLCULO DE CONCRETO CONFINADO (Mander/Kent & Park) ---
def calc_confined_concrete(fc_inicial, epsc0, cobrimento, diam_transversal, diam_longitudinal, npernas_xx, npernas_yy, espacamento_s, b_beam, h_beam, fy_transv, Ast_global):
    """Calcula as propriedades do concreto confinado (fc_c, ec_c, epsu_c)."""
    
    slinha = espacamento_s - diam_transversal
    hx = b_beam - 2 * cobrimento - diam_transversal   
    hy = h_beam - 2 * cobrimento - diam_transversal   
    
    wlinha_x = ((hx - diam_transversal - 2 * diam_longitudinal) - (npernas_yy - 2) * diam_longitudinal) / (npernas_yy - 1)
    wlinha_y = ((hy - diam_transversal - 2 * diam_longitudinal) - (npernas_xx - 2) * diam_longitudinal) / (npernas_xx - 1)

    num_wlinha_x = 2 * (npernas_yy - 1)  
    num_wlinha_y = 2 * (npernas_xx - 1) 
        
    rho_cc = Ast_global / (hx * hy) 
    # Área efetiva de concreto confinado
    Ae = (hx * hy - num_wlinha_x * (wlinha_x**2) / 6 - num_wlinha_y * (wlinha_y**2) / 6) * (1 - 0.5 * slinha / hx) * (1 - 0.5 * slinha / hy)
    Acc = (hx * hy) * (1 - rho_cc) 
    ke = Ae / Acc
        
    # Taxas de armadura transversal
    rho_x = (npernas_xx * math.pi * diam_transversal**2 / 4) / (espacamento_s * hy) 
    rho_y = (npernas_yy * math.pi * diam_transversal**2 / 4) / (espacamento_s * hx) 
        
    # Tensão de confinamento lateral
    fy_transv_N_m2 = fy_transv 
    flx = ke * rho_x * fy_transv_N_m2 
    fly = ke * rho_y * fy_transv_N_m2 
    
    q = np.minimum(flx, fly) / np.maximum(flx, fly) 
    
    # Fatores A e B para Mander/Park
    A_mander = 6.8886 - (0.6069 + 17.275 * q) * np.exp(-4.989 * q)
    B_mander = ((4.5) / (5 / A_mander * (0.9849 - 0.6306 * np.exp(-3.8939 * q)) - 0.1)) - 5
    xlinha = (flx + fly) / (2 * fc_inicial)
    k1 = A_mander * (0.1 + 0.9 / (1 + B_mander * xlinha))
    k2 = 5 * k1  # Fator de aumento de deformação
    
    fator_fc_confinado = 1 + k1 * xlinha
    fator_ec_confinado = 1 + k2 * xlinha
    
    fc_confinado = fc_inicial * fator_fc_confinado
    ec_confinado = epsc0 * fator_ec_confinado 
    
    epsc_u_confinado = 0.008 # Deformação última assumida
    
    if fc_confinado <= fc_inicial:
        fc_confinado = fc_inicial
        ec_confinado = epsc0
        epsc_u_confinado = epsc0
        
    return fc_confinado, ec_confinado, epsc_u_confinado

# --- 1. PROPRIEDADES GEOMÉTRICAS E MATERIAIS ---

# Geometria do pilar
B = 0.55 #m
H = 0.55 #m
L = 1.65 #m
A = B*H 
cobrimento = 0.04 # m
rho_concreto = 2500 # kg/m3

# Aço Longitudinal
matTag_aco = 1
Fy = 511e6 #N/m2 resistência ao escoamento
E0 = 200e9 #N/m2 mod de elasticidade
b = 0.01
diam_longitudinal = 0.020 # m
barras_camadas = [5, 2, 2, 5] 
area_s = np.pi * (diam_longitudinal**2) / 4 
Ast = sum(barras_camadas) * area_s 
ops.uniaxialMaterial('Steel01', matTag_aco, Fy, E0, b)

# Aço Transversal e Parâmetros de Concreto
diam_transversal = 0.012 # m
fy_transv = 325e6 # N/m2 (325 MPa)
npernas_yy = npernas_xx = 4
fc_inicial = 32e6 # N/m2
epsc0_nc = 0.00208

# Parâmetros padrão para Concrete02
rat = 0.6  # Fator de redução de rigidez de recarregamento
ft = 3.0e6 # Resistência à tração (exemplo: 3 MPa)
Ets = 0.0  # Módulo de Elasticidade após a tração

# --- Concreto Não Confinado (Concrete02) ---
matTag_NC = 2
fpc_nc = 32e6      
epsc0_nc = 0.00208 
fpcu_nc = 0.2 * fpc_nc 
epsU_nc = 0.0045 
# Sintaxe CORRIGIDA: fpc e epsc0 POSITIVOS; fpcu e epscu NEGATIVOS
ops.uniaxialMaterial('Concrete02', matTag_NC, fpc_nc, epsc0_nc, -fpcu_nc, -epsU_nc, rat, ft, Ets) 

# --- Concreto Confinado 1 (S1 = 0.11m) - Região Crítica ---
espacamento_s_1 = 0.11 # m
fc_c1, ec_c1, epsu_c1 = calc_confined_concrete(fc_inicial, epsc0_nc, cobrimento, diam_transversal, diam_longitudinal, npernas_xx, npernas_yy, espacamento_s_1, B, H, fy_transv, Ast)
matTag_C1 = 3
ops.uniaxialMaterial('Concrete02', matTag_C1, fc_c1, ec_c1, -fc_c1 * 0.20, -epsu_c1, rat, ft, Ets)

# --- Concreto Confinado 2 (S2 = 0.22m) - Região Meio/Topo ---
espacamento_s_2 = 0.22 # m
fc_c2, ec_c2, epsu_c2 = calc_confined_concrete(fc_inicial, epsc0_nc, cobrimento, diam_transversal, diam_longitudinal, npernas_xx, npernas_yy, espacamento_s_2, B, H, fy_transv, Ast)
matTag_C2 = 4
ops.uniaxialMaterial('Concrete02', matTag_C2, fc_c2, ec_c2, -fc_c2 * 0.20, -epsu_c2, rat, ft, Ets)


# --- 3. NÓS E RESTRIÇÕES ---
node_base = 1
node_topo = 4
L_seg = L / 3.0
ops.node(node_base, 0, 0) 
ops.node(2, 0, L_seg) 
ops.node(3, 0, 2 * L_seg) 
ops.node(node_topo, 0, L) 

ops.fix(node_base, 1, 1, 1) 

# --- 4. SEÇÃO FIBER ---
secTag = 1
ops.section('Fiber', secTag)
c = cobrimento
# Confinado (Núcleo)
ops.patch('quad', matTag_C1, 20, 20, (-H/2)+c, (-B/2)+c, (H/2)-c, (-B/2)+c, (H/2)-c, (B/2)-c, (-H/2)+c, (B/2)-c) 
# Não Confinado (Cobrimento)
ops.patch('rect', matTag_NC, 20, 1, (-H/2), (B/2)-c, (H/2), (B/2))
ops.patch('rect', matTag_NC, 20, 1, (-H/2), (-B/2), (H/2), (-B/2)+c)
ops.patch('rect', matTag_NC, 1, 20, (-H/2), (-B/2)+c, (-H/2)+c, (B/2)-c)
ops.patch('rect', matTag_NC, 1, 20, (H/2)-c, (-B/2)+c, (H/2), (B/2)-c)
# Aço Longitudinal
ops.layer('straight', matTag_aco, 4, area_s, (H/2)-c, (B/2)-c, (H/2)-c, (-B/2)+c)
ops.layer('straight', matTag_aco, 2, area_s, (H-2*c)/6, (B/2)-c, (H-2*c)/6, (-B/2)+c)
ops.layer('straight', matTag_aco, 2, area_s, -(H-2*c)/6, (B/2)-c, -(H-2*c)/6, (-B/2)+c)
ops.layer('straight', matTag_aco, 4, area_s, (-H/2)+c, (B/2)-c, (-H/2)+c, (-B/2)+c)


# --- 5. ELEMENTOS E INTEGRAÇÃO ---
transfTag = 1
ops.geomTransf('PDelta', transfTag) 
integrationTag = 1 
ops.beamIntegration('Lobatto', integrationTag, secTag, 6)

# Elementos 
ops.element('dispBeamColumn', 1, node_base, 2, transfTag, integrationTag) 
ops.element('dispBeamColumn', 2, 2, 3, transfTag, integrationTag) 
ops.element('dispBeamColumn', 3, 3, node_topo, transfTag, integrationTag) 

# --- 6. RECORDER (FORÇA x DESLOCAMENTO) ---
dof_lateral = 1
ops.recorder('Node', '-file', 'Col_Disp.out', '-time', '-node', node_topo, '-dof', dof_lateral, 'disp')
ops.recorder('Node', '-file', 'Col_Forca.out', '-time', '-node', node_base, '-dof', dof_lateral, 'reaction') 

# --- 7. CARGA AXIAL E PESO PRÓPRIO ---
ops.timeSeries("Constant", 1) 
ops.pattern("Plain", 1, 1) 

# Carga axial (0.1 * f'c * A)
P_axial_calc = -(0.1 * fc_inicial * A) 
ops.load(node_topo, 0.0, P_axial_calc, 0.0) 

# Peso próprio distribuído
peso_proprio_por_no = -rho_concreto * A * 9.81 * L_seg
ops.load(2, 0.0, peso_proprio_por_no, 0.0) 
ops.load(3, 0.0, peso_proprio_por_no, 0.0) 
ops.load(node_topo, 0.0, peso_proprio_por_no / 2.0, 0.0) 

# Análise Estática para Carga Axial
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Transformation')
ops.algorithm('Newton')
ops.test('NormUnbalance', 1e-05, 20, 0) # Tolerância 1e-05, 20 iterações
num_passos_axial = 10
ops.integrator('LoadControl', 1.0 / num_passos_axial) # Aplica 1/10 da carga por passo
ops.analysis('Static')
ok = ops.analyze(1)
if ok != 0:
    print("⚠️ Falha ao aplicar a Carga Axial.")
    exit()
ops.loadConst('-time', 0.0)
print("✅ Carga Axial e Peso Próprio aplicados e congelados.")


# --- 8. ANÁLISE CÍCLICA (DISPLACEMENTCONTROL) - UNIDADE 5 ---
Delta_y_M = 0.0123 
node_topo = 4 
dof_lateral = 1    

Mu_picos = [0.75, -0.75, 2, -2, 2, -2, 4, -4, 4, -4, 6, -6, 6, -6]
Disp_Alvo = [0.0] 
for mu in Mu_picos:
    Disp_Alvo.append(mu * Delta_y_M)

dD = 0.00005 # Incremento de deslocamento (0.05 mm)

# Define a 'dummy load' para o DisplacementControl atuar (Pattern 2)
ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
ops.load(node_topo, 1.0, 0.0, 0.0) 

# Configurações de Análise Estática para Carregamento Cíclico
ops.analysis('Static')
# *** AJUSTE CHAVE: USANDO ALGORITMO NEWTONLINE SEARCH PARA MAIOR ROBUSTEZ ***
ops.algorithm('NewtonLineSearch') 
ops.test('NormUnbalance', 1e-04, 10, 0) # Tolerância relaxada

# Loop principal (DisplacementControl)
for i in range(len(Disp_Alvo) - 1):
    D_start = Disp_Alvo[i]
    D_end = Disp_Alvo[i+1]
    dDelta = D_end - D_start
    
    if abs(dDelta) > 1e-6:
        nSteps = int(abs(dDelta) / dD) 
        if nSteps == 0: nSteps = 1
            
        dU_increment = dDelta / nSteps
        ops.integrator('DisplacementControl', node_topo, dof_lateral, dU_increment)
        
        print(f"Segmento {i+1}: {D_start:.4f} m -> {D_end:.4f} m em {nSteps} passos.")
        
        ok = ops.analyze(nSteps) 
        
        # Recuperação de Convergência (Tenta BFGS, se NewtonLineSearch falhar)
        if ok != 0:
            print(f"*** ⚠️ Falha na convergência no segmento {i+1} com NewtonLineSearch. Tentando BFGS... ***")
            ops.algorithm('BFGS')
            ops.integrator('DisplacementControl', node_topo, dof_lateral, dU_increment / 10.0) 
            
            # Tenta analisar o restante dos passos
            nSteps_analisados = ops.analysis('getStatus', 'numSteps') 
            nSteps_restante = nSteps - nSteps_analisados
            
            if nSteps_restante > 0:
                ok = ops.analyze(nSteps_restante)
            
            if ok != 0:
                print(f"--- ❌ Falha dupla no passo {i+1}. Interrompendo. ---")
                break
            # Retorna ao algoritmo principal após a recuperação
            ops.algorithm('NewtonLineSearch') 
            
print("✅ Análise Cíclica Concluída.")


# --- 9. PÓS-PROCESSAMENTO E PLOTAGEM ---
try:
    # Carrega dados ignorando a coluna de tempo (usecols=1)
    disp_data = np.loadtxt('Col_Disp.out', usecols=1)
    force_data = np.loadtxt('Col_Forca.out', usecols=1)
    
    # Plota a curva Histerética
    plt.figure(figsize=(10, 6))
    plt.plot(disp_data, - force_data / 1000, 'b-', linewidth=1) # Dividido por 1000 para Força em kN
    plt.xlabel('Deslocamento do Topo (m)')
    plt.ylabel('Força Lateral (kN)')
    plt.title('Curva Histerética (Força x Deslocamento) - Unidade 5')
    plt.grid(True)
    plt.show()
    
    # Salva os dados
    results_df = pd.DataFrame({
        'Deslocamento (m)': disp_data,
        'Força Lateral (N)': force_data
    })
    results_df.to_excel('resultados_histerese_unidade5.xlsx', index=False, engine='openpyxl')
    print("✅ Arquivo 'resultados_histerese_unidade5.xlsx' criado com sucesso!")
    
except FileNotFoundError:
    print("\n❌ Arquivos de Recorder não encontrados. Verifique a saída do console para erros de convergência.")
except Exception as e:
    print(f"\n❌ Erro durante o pós-processamento: {e}")