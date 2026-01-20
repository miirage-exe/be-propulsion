import matplotlib.pyplot as plt
import numpy as np

# Données de vol
M0 = 0.8
altitude = 36000/3.281
# Ps0 = 227* 10*10
# Ts0 = 217

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


# def motor(fluidStateInit):
#   global alpha
#   print(fluidStateInit)
#   state = tuyereIn(fluidStateInit) # Tuyere 1
#   state_2 = state
#   state = compressor(OPRfan, rendementFan, state) # Fan
#   stateFan = state
#   stateSortieFluxSecondaire = tuyereSortie(stateFan)

#   state = compressor(OPRbp, rendementComp, state) # Compresseur BP
#   stateCompBP = state
#   state = compressor(OPRhp, rendementComp, state) # Compresseur HP
#   stateCompHP = state
#   state = combustionChamber(state, TtSortieCC)
#   state = turbineHP(rendementTurbHP, stateCompBP[1], stateCompHP[1], state)
#   state = turbineBP(rendementTurbBP, stateFan[1], stateCompBP[1], state_2[1], state)
#   stateSortieFluxPrincipal = tuyereSortie(state)
#   print("Sortie de la tuyere: ",state)

#   # Calculs

#   mach_19 = getMach(stateSortieFluxSecondaire[0], Ps_0, gamma)
#   mach_9 = getMach(stateSortieFluxPrincipal[0], Ps_0, gamma_postcomb)
#   mach_0 = getMach(fluidStateInit[0], Ps_0, gamma)

#   vitesse_19 = mach_19*vitesseSon(mach_19, stateSortieFluxSecondaire[1], gamma, r)
#   vitesse_0 = mach_0*vitesseSon(mach_0, fluidStateInit[1], gamma, r)
#   vitesse_9 = mach_9*vitesseSon(mach_9, stateSortieFluxPrincipal[1], gamma_postcomb, r_postcomb)
#   print(vitesse_19, vitesse_9)


#   debit_principal=throttleObj/((1+alpha)*vitesse_9 + lmbda*vitesse_19 - (1+lmbda)*vitesse_0)
#   surface_fan = getSurface(0.6, r, gamma, debit_principal*(1+lmbda), state_2)
#   print(surface_fan)
#   r_min = (debit_principal*(1+lmbda)/3.141592*0.6)**(1/2)
# Main ---------------------------------

# g
# --- Nouvelle fonction de calcul du rendement ---
def calculer_rendement(current_opr, current_Tt4, stateInit, Ps_ambiant):
    global alpha
    current_OPRbp = current_opr / (OPRhp * OPRfan)
    state = tuyereIn(stateInit) 
    state_2 = state
    state = compressor(OPRfan, rendementFan, state) # Fan
    stateFan = state
    stateSortieFluxSecondaire = tuyereSortie(stateFan)
    state = compressor(current_OPRbp, rendementComp, state) # Compresseur BP
    stateCompBP = state
    state = compressor(OPRhp, rendementComp, state) # Compresseur HP
    stateCompHP = state
    state = combustionChamber(state, current_Tt4)
    
    state = turbineHP(rendementTurbHP, stateCompBP[1], stateCompHP[1], state)
    state = turbineBP(rendementTurbBP, stateFan[1], stateCompBP[1], state_2[1], state)
    stateSortieFluxPrincipal = tuyereSortie(state)

    # 3. Calculs des vitesses
    mach_19 = getMach(stateSortieFluxSecondaire[0], Ps_ambiant, gamma)
    mach_9 = getMach(stateSortieFluxPrincipal[0], Ps_ambiant, gamma_postcomb)
    mach_0 = getMach(stateInit[0], Ps_ambiant, gamma)

    vitesse_19 = mach_19 * vitesseSon(mach_19, stateSortieFluxSecondaire[1], gamma, r)
    vitesse_0 = mach_0 * vitesseSon(mach_0, stateInit[1], gamma, r)
    vitesse_9 = mach_9 * vitesseSon(mach_9, stateSortieFluxPrincipal[1], gamma_postcomb, r_postcomb)

    # --- CALCULS DES RENDEMENTS ---
    
    # 1. Puissance Chimique (Input)
    puissance_chimique = alpha * Pk
    
    # 2. Poussée Spécifique (Fsp)
    # F = (m_prim + m_fuel)*V9 + m_sec*V19 - m_total*V0
    # On travaille par kg d'air entrant, donc m_total = 1+lambda
    Fsp = (1 + alpha) * vitesse_9 + lmbda * vitesse_19 - (1 + lmbda) * vitesse_0
    
    # 3. Puissance Propulsive (Utile)
    puissance_propulsive = Fsp * vitesse_0
    
    if puissance_chimique == 0: return 0, 0, 0

    # Calcul des 3 rendements
    eta_global = puissance_propulsive / puissance_chimique
    
    # Note: On peut aussi calculer le thermique et le propulsif séparément
    delta_ec = 0.5 * ((1 + alpha) * (vitesse_9**2) + lmbda * (vitesse_19**2) - (1 + lmbda) * (vitesse_0**2))
    eta_thermique = delta_ec / puissance_chimique
    eta_propulsif = puissance_propulsive / delta_ec if delta_ec != 0 else 0
    
    return eta_global, eta_thermique, eta_propulsif, (Fsp/12)


plage_alt = np.linspace(30_000, 38_000, 40) 


# motor(stateInit)
# --- 1. Modèle Atmosphère Standard (ISA) ---
def get_isa_conditions(altitude_ft):
    # Conversion ft -> mètres
    h = altitude_ft * 0.3048
    
    # Niveau de la mer
    T0_isa = 288.15
    P0_isa = 101325
    
    # Gradient thermique (Troposphère jusqu'à 11km)
    
    Ts = T0_isa - 0.0065 * h
    Ps = P0_isa * (Ts / T0_isa)**5.2561
    return Ps, Ts


# plt.figure(figsize=(10, 6))

rendements_globaux = []
rendements_th = []
rendements_prop = []
force = []
for alt in plage_alt:

    Ps_0, Ts_0 = get_isa_conditions(alt)

    print(Ps_0)

    Pt0 = Ps_0 * (1+ (gamma-1)/2*(M0*M0))**(gamma/(gamma-1))
    Tt0 = Ts_0 * (1+ (gamma-1)/2*(M0*M0))
    stateInit = [Pt0, Tt0, gamma, r, cp]

    eta_g, eta_th, eta_prop, fesp = calculer_rendement(40, 1600, stateInit, Ps_0)
    rendements_globaux.append(eta_g)
    rendements_th.append(eta_th)
    rendements_prop.append(eta_prop)
    force.append(fesp)
    # for crtlmbda in plage_lambda:
    #     # Mise à jour des variables globales pour la simulation
    #     OPRfan = taux_fan
    #     lmbda = crtlmbda
        
        # Appel de la fonction (OPR global fixé à 40, Tt4 à 1600)
    
# plt.plot(plage_alt, rendements_globaux, label=f'Rendement global')
# plt.plot(plage_alt, rendements_th, label=f'Rendement chimique')
# plt.plot(plage_alt, rendements_prop, label=f'Rendement propulsion')

plt.plot(plage_alt, force)

plt.title("Poussée spécifique en fonction de l'altitude")
plt.xlabel("Altitude (ft)")
plt.ylabel("Poussée spécifique (N/(kg/s))")
plt.grid(True)
# plt.legend()
plt.show()