import math

injection_liquid = "Freon-114-B2"

F_1 = 50.0e3
F_2 = 24.0e3

I_sp = 215.0
k = 1.2
g0 = 9.81

D_throat_1 = 86.0e-3
D_throat_2 = 84.0e-3

D_exit_1 = 226.0e-3
D_exit_2 = 247.0e-3

massflow_1 = F_1/(I_sp*g0)
massflow_2 = F_2/(I_sp*g0)

P_max = 5.0e6
T_max = 2720.0
alpha_nozzle = 15.0/(180.0/math.pi)
alpha_injector = 60.0/(180.0/math.pi)

alpha_thrust_vector = 10.0/(180.0/math.pi)

F_1_side = F_1*math.sin(alpha_thrust_vector)
F_2_side = F_2*math.sin(alpha_thrust_vector)

R_gas_universal = 8.31445 #J/k.mol
M_gas = (0.397*28.01 + 0.124*44.01 + 0.115*2.01 + 0.238*18.01 + 0.124*28.0134 + 0.1/100.0*17.001 + 0.2/100.0*1.007)/1000.0
R_gas_specific = R_gas_universal/M_gas

V_throat = (k*g0*R_gas_specific*T_max)**0.5
P_t = P_max/((k+1.0)/2.0)**(k/(k-1.0))

x_nozzle_1 = 0.3 #Location of the injection port in the nozzle
x_nozzle_2 = 0.3 #Location of the injection port in the nozzle

alpha_discharge_1 = (180.0 + 25.0)/(180.0/math.pi)
alpha_discharge_2 = (180.0 + 25.0)/(180.0/math.pi)

dt_delay = 0.4 #seconds
t_injection = 8.0 #seconds

L_nozzle_1 = (((D_exit_1 - D_throat_1) / 2.0)/math.sin(alpha_nozzle)) * math.cos(alpha_nozzle)
L_nozzle_2 = (((D_exit_2 - D_throat_2) / 2.0)/math.sin(alpha_nozzle)) * math.cos(alpha_nozzle)

L_injector_1 = L_nozzle_1 * x_nozzle_1
L_injector_2 = L_nozzle_2 * x_nozzle_2

D_x_1 = (D_throat_1/2.0 + (L_injector_1/math.cos(alpha_nozzle)) * math.sin(alpha_nozzle)) * 2.0
D_x_2 = (D_throat_2/2.0 + (L_injector_2/math.cos(alpha_nozzle)) * math.sin(alpha_nozzle)) * 2.0

A_x_1 = math.pi * (D_x_1/2.0)**2.0
A_x_2 = math.pi * (D_x_2/2.0)**2.0

A_t_1 = math.pi * (D_throat_1/2.0)**2.0
A_t_2 = math.pi * (D_throat_2/2.0)**2.0

P_x_1 = P_t
P_x_2 = P_t

Q_1 = 2.0
Q_2 = 2.0

dP = 1.0e1
P_atmosphere = 101.25e3

while Q_1 > 1.0:

    Q_1 = (V_throat * ((k+1.0)/(k-1.0) * (1.0 - (2.0/(k+1.0)) * (P_x_1/P_t) ** ((k-1.0)/k)))) / (V_throat * (P_t/P_x_1)**(1/k) * (A_t_1/A_x_1))

    P_x_1 = P_x_1 - dP

while Q_2 > 1.0:

    Q_2 = (V_throat * ((k+1.0)/(k-1.0) * (1.0 - (2.0/(k+1.0)) * (P_x_2/P_t) ** ((k-1.0)/k)))) / (V_throat * (P_t/P_x_2)**(1/k) * (A_t_2/A_x_2))

    P_x_2 = P_x_2 - dP

print("Total nozzle 1 length: {0} mm".format(L_nozzle_1 * 1000.0))
print("Total nozzle 2 length: {0} mm\n".format(L_nozzle_2 * 1000.0))

print("Injector 1 distance: {0} mm".format(L_injector_1 * 1000.0))
print("Injector 2 distance: {0} mm\n".format(L_injector_2 * 1000.0))

print("Pressure at injection point 1: {0} Pa".format(P_x_1 + P_atmosphere))
print("Pressure at injection point 2: {0} Pa".format(P_x_2 + P_atmosphere))

