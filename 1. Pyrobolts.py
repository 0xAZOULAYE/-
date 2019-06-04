import math

#Pyrotechnic bolt design

mat_name = "304"
mat_E = 200.0*10**9.0
mat_sigma_UTS = 515.0*10**6 #Pa
mat_sigma_yield = 205.0*10**6 #Pa
mat_sigma_fatigue = 140.0*10**6.0 #Pa

mat_sigma_K = 0.0

D_internal = 10.0*10**-3
L_internal = 50.0*10**-3
t_wall = 1.5*10**-3

T_g_petn = 4230.0 + 273.0 #K
R_universal = 8.314 #J/mol.K
M_petn = 0.316 #kg/mol
rho_petn = 1770.0 #kg/m3

E_min_spark = 10*10**-3 #Electric spark energy needed for direct initiation
E_max_spark = 60*10**-3 #Electric spark energy needed for direct initiation

V_battery = 12.0 #V

F_pull = 0.0

#Determine the needed breaking pressure
sigma_hoop = 0.0
sigma_z = 0.0
pressure = 101250.0

safety_factor = 2.0

while sigma_hoop < mat_sigma_UTS and sigma_z < mat_sigma_UTS:

    sigma_hoop = pressure*(D_internal/2.0)/t_wall
    sigma_z = pressure*(D_internal/2.0)/(2.0*t_wall)
    
    pressure = pressure + 1000.0

print("The breaking pressure is: {0} MPa".format(pressure/1000000.0))

#Define the maximum tensile force
A_cross = math.pi*(D_internal/2.0 + t_wall)**2.0 - math.pi*(D_internal/2.0)**2.0
F_pull = A_cross * mat_sigma_fatigue

print("The dynamic tensile loading is: {0} N".format(F_pull))

#Define the needed PETN loading
V_chamber = math.pi*(D_internal/2.0)**2.0*L_internal
m_charge = (pressure*V_chamber)/((R_universal/M_petn)*T_g_petn)

V_charge = m_charge/rho_petn
L_charge = V_charge/(math.pi*(D_internal/2.0)**2.0)

print("The PETN charge loading is: {0} grams".format(m_charge*1000.0))
print("The PETN charge loading is: {0} mm".format(L_charge*1000.0))


#Define the needed primer capacitor
C_detonator = (E_max_spark*2.0)/(V_battery**2.0)
print("The primer capacitor charge loading is: {0} micro-farad".format(C_detonator*1000000.0))


