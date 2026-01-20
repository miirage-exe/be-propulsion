import matplotlib.pyplot as plt
import numpy as np

# Données de vol
M0 = 0.8
altitude = 36000/3.281
Ps0 = 227* 10*10
Ts0 = 217

# Données fluide
gamma = 1.4
gamma_postcomb = 1.33
r = 287
r_postcomb = 291.6
cp = gamma * r /(gamma-1)
cp_postcomb = gamma_postcomb * r_postcomb /(gamma_postcomb-1)
Pk = 42_800_000 # Pouvoir calorifique du kérosène

# Données objectif
throttleObj = 21000
consumptionObj = 0

# Données structurelles
compressionRatio = 1.45
BPR = 11 # Taux de dilution
OPR = 40# Taux de compression global (entre 2 et 3)
OPRhp = 22 # Taux de compression compresseur Haute Pression ( entre 2_5 et 3)
OPRfan = 1.45 #Taux de compression entre 2 et 2_1
OPRbp = OPR/(OPRhp*OPRfan) #Taux de compression entre 2_1 et 2_5
tsetaE = 0.98
TtSortieCC = 1600
lmbda = 11
rRatio = 0.6

rendementComp = 0.9
rendementFan = 0.92
rendementCombustion = 0.99
tsetaCC = 0.95
rendementMeca = 0.98
rendementTurbHP = 0.89
rendementTurbBP = 0.90
tsetaTuyere = 0.98

fluideState = [] # Pt Tt gamma r cp debit
alpha = 0

# Entre 0 et 2
def tuyereIn(fluideState):
  newPt = tsetaE * fluideState[0]
  return [newPt] + fluideState[1:]

def compressor(tauxComp, rendement, fluideState):
  newPt = tauxComp * fluideState[0]
  newTt = (tauxComp)**((fluideState[2]-1)/(fluideState[2]*rendement))*fluideState[1]
  return [newPt, newTt] + fluideState[2:]

def combustionChamber(fluideState, TtTarget):
  global alpha
  newPt = fluideState[0]*tsetaCC
  alpha = (cp_postcomb*TtTarget - fluideState[4]*fluideState[1])/(rendementCombustion*Pk - TtTarget*cp_postcomb)
  newTt = TtTarget
  return [newPt, newTt, gamma_postcomb, r_postcomb, cp_postcomb]

def turbineHP(rendement, TinComp, ToutComp, fluideState):
  newTt = fluideState[1]-(cp/(cp_postcomb*rendementMeca*(1+alpha))*(ToutComp-TinComp))
  pic = (fluideState[1]/newTt)**(fluideState[2]/(rendement*(fluideState[2]-1)))
  newPt = fluideState[0]/pic
  return [newPt, newTt] + fluideState[2:]

def turbineBP(rendement, TinComp, ToutComp, TinFan, fluideState):
  newTt = fluideState[1]-(cp/(cp_postcomb*rendementMeca*(1+alpha))*(ToutComp+lmbda*TinComp -TinFan*(1+lmbda)))
  pic = (fluideState[1]/newTt)**(fluideState[2]/(rendement*(fluideState[2]-1)))
  newPt = fluideState[0]/pic
  return [newPt, newTt] + fluideState[2:]


def tuyereSortie(fluideState):
  newPt = tsetaTuyere*fluideState[0]
  newTt = fluideState[1]
  return [newPt, newTt] + fluideState[2:]

def getMach(pt, ps, gamma):
  return (2/(gamma-1)*((pt/ps)**((gamma-1)/gamma)-1))**(1/2)

def vitesseSon(mach, tt, gamma, r):
  tempStatique=tt/(1+(gamma-1)/2*mach*mach)
  return (gamma*r*tempStatique)**(1/2)

def sigma(mach, r, gamma):
  return (gamma/r)**(1/2)*mach*(1+(gamma-1)/2*mach*mach)**(-(gamma+1)/(2*(gamma-1)))

def getSurface(mach, r, gamma, debit, fluideState):
  sig = sigma(mach, r, gamma)
  return debit/sig*(fluideState[1])**(1/2)/fluideState[0]


def motor(fluidStateInit):
  global alpha
  print(fluidStateInit)
  state = tuyereIn(fluidStateInit) # Tuyere 1
  state_2 = state
  state = compressor(OPRfan, rendementFan, state) # Fan
  stateFan = state
  stateSortieFluxSecondaire = tuyereSortie(stateFan)

  state = compressor(OPRbp, rendementComp, state) # Compresseur BP
  stateCompBP = state
  state = compressor(OPRhp, rendementComp, state) # Compresseur HP
  stateCompHP = state
  state = combustionChamber(state, TtSortieCC)
  state = turbineHP(rendementTurbHP, stateCompBP[1], stateCompHP[1], state)
  state = turbineBP(rendementTurbBP, stateFan[1], stateCompBP[1], state_2[1], state)
  stateSortieFluxPrincipal = tuyereSortie(state)
  print("Sortie de la tuyere: ",state)

  # Calculs

  mach_19 = getMach(stateSortieFluxSecondaire[0], Ps0, gamma)
  mach_9 = getMach(stateSortieFluxPrincipal[0], Ps0, gamma_postcomb)
  mach_0 = getMach(fluidStateInit[0], Ps0, gamma)

  vitesse_19 = mach_19*vitesseSon(mach_19, stateSortieFluxSecondaire[1], gamma, r)
  vitesse_0 = mach_0*vitesseSon(mach_0, fluidStateInit[1], gamma, r)
  vitesse_9 = mach_9*vitesseSon(mach_9, stateSortieFluxPrincipal[1], gamma_postcomb, r_postcomb)
  print(vitesse_19, vitesse_9)


  debit_principal=throttleObj/((1+alpha)*vitesse_9 + lmbda*vitesse_19 - (1+lmbda)*vitesse_0)
  surface_fan = getSurface(0.6, r, gamma, debit_principal*(1+lmbda), state_2)
  print(surface_fan)
  r_min = (debit_principal*(1+lmbda)/3.141592*0.6)**(1/2)
# Main ---------------------------------

plage_lambda = np.linspace(8, 15, 40) # OPR de 20 à 60
plage_opr = [1.1, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40] # Différentes TtSortieCC en Kelvin
# plt.figure(figsize=(10, 6))

Pt0 = Ps0 * (1+ (gamma-1)/2*(M0*M0))**(gamma/(gamma-1))
Tt0 = Ts0 * (1+ (gamma-1)/2*(M0*M0))
stateInit = [Pt0, Tt0, gamma, r, cp]

# motor(stateInit)

# g
# --- Nouvelle fonction de calcul du rendement ---
def calculer_rendement(current_opr, current_Tt4, stateInit):
    global alpha
    
    # 1. Recalculer OPRbp en fonction du OPR global demandé
    current_OPRbp = current_opr / (OPRhp * OPRfan)
    
    # 2. Exécution du cycle
    state = tuyereIn(stateInit) 
    state_2 = state
    
    state = compressor(OPRfan, rendementFan, state) # Fan
    stateFan = state
    stateSortieFluxSecondaire = tuyereSortie(stateFan)

    state = compressor(current_OPRbp, rendementComp, state) # Compresseur BP
    stateCompBP = state
    
    state = compressor(OPRhp, rendementComp, state) # Compresseur HP
    stateCompHP = state
    
    # Passage de la Tt4 cible à la chambre de combustion
    state = combustionChamber(state, current_Tt4)
    
    state = turbineHP(rendementTurbHP, stateCompBP[1], stateCompHP[1], state)
    state = turbineBP(rendementTurbBP, stateFan[1], stateCompBP[1], state_2[1], state)
    stateSortieFluxPrincipal = tuyereSortie(state)
    print("Sortie de la tuyere: ",stateSortieFluxPrincipal)

    # 3. Calculs des vitesses
    mach_19 = getMach(stateSortieFluxSecondaire[0], Ps0, gamma)
    mach_9 = getMach(stateSortieFluxPrincipal[0], Ps0, gamma_postcomb)
    mach_0 = getMach(stateInit[0], Ps0, gamma)

    vitesse_19 = mach_19 * vitesseSon(mach_19, stateSortieFluxSecondaire[1], gamma, r)
    vitesse_0 = mach_0 * vitesseSon(mach_0, stateInit[1], gamma, r)
    vitesse_9 = mach_9 * vitesseSon(mach_9, stateSortieFluxPrincipal[1], gamma_postcomb, r_postcomb)

    # --- CALCULS DES RENDEMENTS ---
    
    # 1. Puissance Chimique (Input)
    puissance_chimique = alpha * Pk
    
    # 2. Poussée (F)
    # F = (m_prim + m_fuel)*V9 + m_sec*V19 - m_total*V0
    F = (1 + alpha) * vitesse_9 + lmbda * vitesse_19 - (1 + lmbda) * vitesse_0

    # 3. Puissance Propulsive (Utile)
    puissance_propulsive = F * vitesse_0
    print("puissance prop:", puissance_propulsive)
    
    if puissance_chimique == 0: return 0, 0, 0

    # Calcul des 3 rendements
    eta_global = puissance_propulsive / puissance_chimique
    
    # Note: On peut aussi calculer le thermique et le propulsif séparément
    delta_ec = 0.5 * ((1 + alpha) * (vitesse_9**2) + lmbda * (vitesse_19**2) - (1 + lmbda) * (vitesse_0**2))
    eta_thermique = delta_ec / puissance_chimique
    eta_propulsif = puissance_propulsive / delta_ec if delta_ec != 0 else 0
    
    return eta_global, eta_thermique, eta_propulsif


# plage_lambda = np.linspace(2, 40, 100) 
# plage_opr_fan = [1.25, 1.35, 1.45, 1.65, 1.85] 

# plt.figure(figsize=(10, 6))

eta_g, eta_th, eta_p = calculer_rendement(40, 1600, stateInit)
    
# plt.plot(plage_lambda, rendements_globaux, label=f'Taux compression fan = {taux_fan}')

# plt.title("Rendement Global en fonction du BPR (Lambda) pour différents taux de compression du fan")
# plt.xlabel("Taux de dilution (BPR/Lambda)")
# plt.ylabel("Rendement Global")
# plt.grid(True)
# plt.legend()
# plt.show()