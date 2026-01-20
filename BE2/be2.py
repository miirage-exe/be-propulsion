import math

T0 = 217
P0 = 22700
M0 = 0.82
epsilon = 0.98
gamma = 1.4
M12 = 0.6
r = 287.1
debit12 = 170
ratio = 0.3
pi_f = 1.45
rendement_f = 0.92

Pt0 = P0 * (1+(gamma-1)/2*M0*M0)**(gamma/(gamma-1))
Tt0 = T0 * (1+(gamma-1)/2*M0*M0)

Pt12 = epsilon * Pt0
Tt12 = Tt0

Ps12 = Pt12/(1+(gamma-1)/2*M12*M12)**(gamma/(gamma-1))
Ts12 = Tt12/(1+(gamma-1)/2*M12*M12)

rho = Ps12/(r*Ts12)

a12 = (gamma*r*Ts12)**(1/2)
v12 = a12*M12

surface = debit12 / (rho*v12)
rmax = (surface / (3.1415 * (1-(ratio)**2)))**(1/2)
rmin = 0.3*rmax

print("Rmax: ", rmax)
print("Rmin: ", rmax*0.3)
print("Dmax: ", 2*rmax)

#Q2
Vw_max = 1.5*a12
Umax = (Vw_max*Vw_max-v12*v12)**(1/2)

thetar_max = Umax/rmax

print("vr_max (rad/s): ", thetar_max)

#Q3
U_pied = thetar_max*rmin
U_tete = thetar_max*rmax
vrel_pied = [v12 , thetar_max*rmin]
vrel_tete = [v12 , thetar_max*rmax]

beta12_pied = math.atan(vrel_pied[1]/vrel_pied[0])
beta12_tete = math.atan(vrel_tete[1]/vrel_tete[0])

print("Beta 12 pied (deg): ", beta12_pied*180/3.141592)
print("Beta 12 tete (deg): ", beta12_tete*180/3.141592)
# On s'impose une incidence nulle au niveau du pied de la pale

#Q4
# En compressant on augmente Pt et donc rho (Ps % à Pt).
# On a bien reflechit et faut que ça diminue

#Q5

# Il suffit de calculer la nouvelle vitesse en 13
Tt13 = Tt12* (pi_f)**((gamma-1)/(gamma*rendement_f))
cp = gamma*r/(gamma-1)
print("cp:", cp)

ht = (Tt13-Tt12)*cp
scalaire = ht # Formule d'Euler
v13 = v12

U_13_pied = thetar_max*rmax*0.4

#alpha = 3.141592 - math.acos(scalaire/(U_13_pied*v13))


#print("alpha", alpha*180/3.1415)

# vx12_pied = v12*math.cos(tmp_pied)
# beta13_pied = (U_pied + vx12_pied*math.tan(tmp_pied))/vx12_pied

# print("Beta 13 pied: ", beta13_pied*180/3.1415)

v_theta = scalaire/U_13_pied
print("v_theta", v_theta)
print("U_pied", U_13_pied)
alpha = math.atan(v_theta/v13)
print("alpha", alpha*180/3.1415)
print("v13", v13)
W_theta_13 = v_theta -  U_13_pied
print("W:",W_theta_13)
beta13_pied = math.atan(W_theta_13/v13)
print("Beta  13 (pied) : ",  beta13_pied*180/3.1415)



v_theta = scalaire/U_tete
print("v_theta", v_theta)
print("U_pied", U_tete)
alpha = math.atan(v_theta/v13)
print("alpha", alpha*180/3.1415)
print("v13", v13)
W_theta_13 = v_theta -  U_tete
print("W:",W_theta_13)
beta13_tete = math.atan(W_theta_13/v13)
print("Beta  13 (tete) : ",  beta13_tete*180/3.1415)

# Q6
# Faut que ce soit alpha

# Q7
delta_beta_pied = beta12_pied - beta13_pied
delta_beta_tete = beta12_tete - beta13_tete  

print("Diff. Beta gros orteil (deb): ", delta_beta_pied*180/3.141592)
print("Diff. Beta tete (deb): ", delta_beta_tete*180/3.141592)