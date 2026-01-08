

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

def combustionChamber(fluideState):
  global alpha
  newPt = fluideState[0]*tsetaCC
  alpha = (cp_postcomb*TtSortieCC - fluideState[4]*fluideState[1])/(rendementCombustion*Pk - TtSortieCC*cp_postcomb)
  newTt = TtSortieCC
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

# Main
Pt0 = Ps0 * (1+ (gamma-1)/2*(M0*M0))**(gamma/(gamma-1))
Tt0 = Ts0 * (1+ (gamma-1)/2*(M0*M0))
state = [Pt0, Tt0, gamma, r, cp]
print(state)
state = tuyereIn(state)
stateIn = state
state = compressor(OPRfan, rendementFan, state) # Fan
stateFan = state
state = compressor(OPRbp, rendementComp, state) # Compresseur BP
stateCompBP = state
state = compressor(OPRhp, rendementComp, state) # Compresseur HP
stateCompHP = state
state = combustionChamber(state)
print(alpha)
state = turbineHP(rendementTurbHP, stateCompBP[1], stateCompHP[1], state)
print("Sortie turbine HP: ", state)
state = turbineBP(rendementTurbBP, stateFan[1], stateCompBP[1], stateIn[1], state)
print("Sortie turbine BP: ",state)