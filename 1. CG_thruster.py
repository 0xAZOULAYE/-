import math

#Cold gas thruster design

m_total = 500.0 #kg
D_outer = 500.0*10**-3 #m

omega_max = 1.8 #dps
t_acceleration = 0.1 #seconds

I_airframe = m_total*(D_outer/2)**2
alpha_rotation = omega_max/t_acceleration

T_rotation = I_airframe*alpha_rotation
number_of_thrusters = 2.0
F_total = T_rotation/(D_outer/2.0)
F_thruster = F_total/number_of_thrusters

number_of_firings = 10.0

P_atm = 101.25e3
P_max = P_atm + 1.0e6 #P
P_tank = 35.0e6 #P

V_max = 1.0e-3
R_air = 287.05
T_tank = 293

#Determine the density of the air
rho_air = (P_max*V_max)/(R_air*T_tank)

gas_type = "compressed air"

#Determine the throat velocity
gamma = 1.4
g_c = 9.81
V_throat = (gamma*g_c*R_air*T_tank)**0.5

#Determine the exit velocity
V_exit = ((2*gamma*g_c*R_air*T_tank)/(gamma-1)*(1-(P_atm/P_max)**((gamma-1)/gamma)))**0.5

#Determine the mass flow per thruster
m_flow = F_thruster/V_exit

#Determine the specific impulse
I_isp = V_exit/g_c
I_total = t_acceleration*F_thruster*number_of_firings*number_of_thrusters

#Determine the throat area from the characteristic velocity
C_star = ((gamma*g_c*R_air*T_tank)**0.5)/ (gamma * ((2/gamma+1)**((gamma+1)/(gamma-1))))
w_propellant = m_flow*g_c

A_throat = (C_star*w_propellant)/(P_max*g_c)
D_throat = 2.0*(A_throat/math.pi)**0.5

#Determine the exit area from the thrust coefficient
C_f = F_thruster/(P_max*A_throat)

#A_exit = (C_f - (((2*gamma**2)/(gamma-1))*(2/(gamma+1))**((gamma+1)/(gamma-1))*(1-(P_exit/P_max)**((gamma-1)/gamma)))**0.5)/(((P_exit-P_atm)/(P_max)))*A_throat

P_exit = P_atm/100.0

A_exit = (( ((gamma*(2/(gamma+1)))**((gamma+1)/(gamma-1)))**0.5)/((P_exit/P_max)**(1/gamma) * ((2*gamma)/(gamma-1) * (1-(P_exit/P_max)**((gamma-1)/gamma) ) )**0.5))*A_throat
D_exit = 2.0*(A_exit/math.pi)**0.5

#Determine the total propellant mass
mass_propellant = I_total/I_isp
volume_propellant = (mass_propellant*R_air*T_tank)/P_tank

D_tank_propellant = 250.0e-3
H_tank_propellant = volume_propellant/(math.pi*(D_tank_propellant/2.0)**2.0)

#Nozzle entry angle
alpha_entry = 45.0

#Nozzle exit angle
alpha_exit = 15.0
L_exit_nozzle = ((D_exit/2.0 - D_throat/2.0)/math.sin(alpha_exit/(180.0/math.pi)))*math.cos(alpha_exit/(180.0/math.pi))

#Nozzle throat radius
throat_radius = 1.5*(D_throat/2.0)

print("Thruster force: {0} N".format(F_thruster))
print("Throat velocity: {0} m/s".format(V_throat))
print("Exit velocity: {0} m/s".format(V_exit))
print("Exit mass flow: {0} kg/s\n".format(m_flow))

print("Exit Isp: {0} s".format(I_isp))
print("Total impulse per thruster: {0} Ns".format(I_total))
print("Total propellant mass: {0} kg".format(mass_propellant))

print("Propellant tank diameter: {0} m".format(D_tank_propellant))
print("Propellant tank height: {0} m".format(H_tank_propellant))

print("Total propellant volume: {0} m3\n".format(volume_propellant))


print("Throat diameter: {0} mm".format(D_throat*1000.0))
print("Throat radius: {0} mm".format(throat_radius*1000.0))

print("Exit diameter: {0} mm".format(D_exit*1000.0))
print("Exit length: {0} mm".format(L_exit_nozzle*1000.0))

