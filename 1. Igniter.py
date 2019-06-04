import math

"""
    Primary propellant composition:
    40% Boron
    60% KNO3
"""

Propellant_primary = "Zr-PETN-PbO2-BaNO3"
Propellant_secondary = "NG-NC-DEP-2NDPA"

R_rod = 176.0/2000.0
W_tube = 45.0e-3
R_heater = 5.0e-3
R_igniter = R_rod - 5.0/1000.0

pressure = 3700000.0
R_universal = 8.13
M_bkno3 = 111.9142/1000.0
T_g_bkno3 = 3070.0

L_grain_1 = 6000.0/1000.0 #m
L_grain_2 = 814.0/1000.0 #m
L_grain_3 = 535.0/1000.0 #m

Overfill_factor = 2.0
V_max = (math.pi*(R_rod+W_tube)**2.0 - math.pi*R_rod**2.0) * L_grain_1 #Maximum fill volume
V_propellant = V_max*Overfill_factor

T_ignition_primary = 351.0 + 273.0 #Ignition propellant ignition temperature
T_gas_primary = 3000.0 #K
T_ignition_secondary = 170.0 + 273.0 #Main propellant ignition temperature

igniter_reaction_rate = 0.61*10**-3 #/s

rho_boron = 2370.0 #kg/m3
rho_kno3 = 2110.0 #kg/m3
rho_bkno3 = 0.4*rho_boron + 0.6*rho_kno3


#Define the needed PETN loading
V_chamber = V_propellant
m_charge = (pressure*V_chamber)/((R_universal/M_bkno3)*T_g_bkno3)

V_charge = m_charge/rho_bkno3
L_charge = V_charge/(math.pi*(R_igniter)**2.0 - math.pi*(R_heater)**2.0)

print("Igniter grain mass: {0} kg".format(m_charge))
print("Igniter grain length: {0} mm".format(L_charge*1000.0))
print("Igniter grain diameter: {0} mm".format(2.0*R_igniter*1000.0))
print("Igniter grain heater diameter: {0} mm".format(2.0*R_heater*1000.0))