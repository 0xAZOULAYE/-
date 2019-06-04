"""
//This script is meant to dimension the internal grain and nozzle of the solid rocket motor
//The following aspects are calculated:
    // - Grain core diameter
    // - Grain outer diameter
    // - Grain length
    // - Grain shape
    // - Grain burn profile
    // - Nozzle entry diameter
    // - Nozzle throat diameter
    // - Nozzle exit diameter
    // - Nozzle throat radius
    // - Nozzle length

    	Composition:		1. NC:		53.0	% 
							2. NG:		40.5	%
							3. DEP:		4.0		%
							4. 2NDPA:	2.5		%
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from decimal import *

getcontext().prec = 7

m_payload = 50.0
n_stages = 1.0

frac_payload = 10.0/100.0
frac_structural = 20.0/100.0
frac_propellant = 1.0 - frac_payload - frac_structural


m_total = 620.0 #1/frac_payload * m_payload
m_propellant = 370.0 #m_total * frac_propellant
m_structural = 200.0 #m_total * frac_structural

m_capsule = m_payload * (1 + frac_structural)
m_stages = m_total - m_capsule

m_individual_stage_1 = m_stages - m_capsule
m_individual_stage_structural_1 = m_individual_stage_1

#m_individual_stage_2 = m_stages/16.0*4.0
#m_individual_stage_structural_2 = m_individual_stage_2*frac_structural

#m_individual_stage_3 = m_stages/16.0*1.0
#m_individual_stage_structural_3 = m_individual_stage_3*frac_structural

Cd_z = 0.03
Cd_r = 1.0

D_body = 442.0/1000.0
D_capsule = 442.0/1000.0


def main():

    rho_NC = 770.0
    rho_NG = 1600.0
    rho_DEP = 1120.0
    rho_2NDPA = 1300.0

    rho_grain_NG_NC = rho_NC*0.53 + rho_NG*0.405 + rho_DEP*0.04 + rho_2NDPA*0.025

    print("Rho grain: {0}".format(rho_grain_NG_NC))

    T_flame = 2720.0                #K

    I_sp = 215.0
    g0 = 9.81

    F_thrust_1 = 50000.0             #N
    #F_thrust_2 = 22000.0             #N
    #F_thrust_3 = 17000.0             #N

    t_thrust = 24.0                 #seconds

    D_max = 452.0/1000.0             #m
    L_max = 4.0                     #m
    m_max = 397.0                   #kg

    motor_OD = D_max
    motor_ID = 40.0/1000.0

    #Determine the motor operation conditions
    P_max = 4000000.0

    P_max_1 = 3700000.0
    #P_max_2 = 2500000.0
    #P_max_3 = 1900000.0

    #NG-NC propellant
    n_motor = 0.64
    a_motor = (2.3*10**-3)/((1.0*10.0**6)**n_motor)

    H_exp_min = 3.47*10**6 #J/kg
    H_exp_max = 4.59*10**6 #J/kg
    

    r_motor_1 = burn_rate(P_max_1, a_motor, n_motor)
    #r_motor_2 = burn_rate(P_max_2, a_motor, n_motor)
    #r_motor_3 = burn_rate(P_max_3, a_motor, n_motor)


    #Grain geometry = rod and tube
    R_rod_support = 5.0e-3
    W_tube = 45.0e-3
    R_rod = 88.0e-3 #Starting diameter rod
    R_tube = R_rod + W_tube #Starting diameter tube
    
    #n_stages, t_stage = number_of_stages(r_motor, R_rod, t_thrust, R_rod_support)

    n_stages, t_stage_1 = number_of_stages(r_motor_1, R_rod, t_thrust, R_rod_support)
    #n_stages, t_stage_2 = number_of_stages(r_motor_2, R_rod, t_thrust, R_rod_support)
    #n_stages, t_stage_3 = number_of_stages(r_motor_3, R_rod, t_thrust, R_rod_support)



    motor_OD_final = (R_tube + r_motor_1*t_stage_1) * 2.0

    m_flow_1 = (F_thrust_1/(I_sp*g0))
    #m_flow_2 = (F_thrust_2/(I_sp*g0))
    #m_flow_3 = (F_thrust_3/(I_sp*g0))

    m_flow_total = 0.0

    #Compute the initial grain surface area
    A_burn_1 = m_flow_1/(rho_grain_NG_NC*r_motor_1)
    #A_burn_2 = m_flow_2/(rho_grain_NG_NC*r_motor_2)
    #A_burn_3 = m_flow_3/(rho_grain_NG_NC*r_motor_3)

    #Determine the length of the grain
    L_grain_1 = A_burn_1/(2.0*math.pi*R_rod + 2.0*math.pi*R_tube)
    #L_grain_2 = A_burn_2/(2.0*math.pi*R_rod + 2.0*math.pi*R_tube)
    #L_grain_3 = A_burn_3/(2.0*math.pi*R_rod + 2.0*math.pi*R_tube)

    erosive_factor_1 = [0.0]
    #erosive_factor_2 = [0.0]
    #erosive_factor_3 = [0.0]

    L_grain_total = L_grain_1 #+ L_grain_2 + L_grain_3
    V_grain_total = (((math.pi * R_rod**2.0) - (math.pi * R_rod_support**2.0)) + ((math.pi * (motor_OD_final/2.0)**2.0) - (math.pi * R_tube**2.0))) * (L_grain_1)
    m_grain_total = V_grain_total*rho_grain_NG_NC

    #Plot the massflow rate over the thrust period
    P_exit_1 = 101250.0
    #P_exit_2 = 80000.0
    #P_exit_3 = 48200.0

    n_steps = 1000.0
    dt = 1.0/n_steps
    t_burn = 0.0
    timestamp_1 = [0.0]
    #timestamp_2 = [0.0]
    #timestamp_3 = [0.0]

    n_counter = 1.0
    
    A_burn_intermediate_1 = A_burn_1
    #A_burn_intermediate_2 = A_burn_2
    #A_burn_intermediate_3 = A_burn_3

    m_flow_intermediate_1 = [0.0]
    #m_flow_intermediate_2 = [0.0]
    #m_flow_intermediate_3 = [0.0]

    F_thrust_intermediate_1 = [0.0]
    #F_thrust_intermediate_2 = [0.0]
    #F_thrust_intermediate_3 = [0.0]

    P_c_intermediate_1 = [0.0]
    #P_c_intermediate_2 = [0.0]
    #P_c_intermediate_3 = [0.0]

    r_motor_intermediate_1 = [0.0]
    #r_motor_intermediate_2 = [0.0]
    #r_motor_intermediate_3 = [0.0]

    M_exhaust_1 = [0.0]
    #M_exhaust_2 = [0.0]
    #M_exhaust_3 = [0.0]

    rod_stage_1 = R_rod
    #rod_stage_2 = R_rod
    #rod_stage_3 = R_rod

    Isp_stage_1 = [0.0]
    #Isp_stage_2 = [0.0]
    #Isp_stage_3 = [0.0]

    R_gas_universal = 8.31445 #J/k.mol
    M_gas = (0.397*28.01 + 0.124*44.01 + 0.115*2.01 + 0.238*18.01 + 0.124*28.0134 + 0.1/100.0*17.001 + 0.2/100.0*1.007)/1000.0 #(R_gas_universal * T_flame * rho_grain_NG_NC)/P_max
    R_gas_specific = R_gas_universal/M_gas

    k = 1.2
    Cp = R_gas_specific/(1.0 - 1.0/k)

    A_throat_1 = m_flow_1/(((((2/(k+1))**((k+1)/(k-1)))**0.5) / ((k*R_gas_specific*T_flame)**0.5))*P_max_1*k)*0.65
    D_throat_1 = ((A_throat_1/math.pi)**0.5)*2.0

    #Define the exit diameter
    A_exit_1 = A_throat_1 * (((k* (2.0/(k+1.0) )**((k+1)/(k-1)))**0.5) / ((P_exit_1/P_max_1)**(1.0/k) * (((2.0*k) / (k-1)) * (1.0 - (P_exit_1/P_max_1) ** ((k-1)/k)) )**0.5))*4.0
    D_exit_1 = ((A_exit_1/math.pi)**0.5)*2.0

    print("Motor 1 throat diameter: {0} mm".format(D_throat_1*1000.0))
    print("Motor 1 exit diameter: {0} mm".format(D_exit_1*1000.0))
    
    """
    A_throat_2 = (m_flow_2/(((((2/(k+1))**((k+1)/(k-1)))**0.5) / ((k*R_gas_specific*T_flame)**0.5))*P_max_2*k))
    D_throat_2 = ((A_throat_2/math.pi)**0.5)*2.0

    #Define the exit diameter
    A_exit_2 = A_exit_1 #A_throat_2 * (((k* (2.0/(k+1.0) )**((k+1)/(k-1)))**0.5) / ((P_exit_2/P_max_2)**(1.0/k) * (((2.0*k) / (k-1)) * (1.0 - (P_exit_2/P_max_2) ** ((k-1)/k)) )**0.5))
    D_exit_2 = ((A_exit_2/math.pi)**0.5)*2.0

    print("Motor 2 throat diameter: {0} mm".format(D_throat_2*1000.0))
    print("Motor 2 exit diameter: {0} mm".format(D_exit_2*1000.0))

    A_throat_3 = (m_flow_3/(((((2/(k+1))**((k+1)/(k-1)))**0.5) / ((k*R_gas_specific*T_flame)**0.5))*P_max_3*k))
    D_throat_3 = ((A_throat_3/math.pi)**0.5)*2.0

    #Define the exit diameter
    A_exit_3 = A_exit_1 #A_throat_3 * (((k* (2.0/(k+1.0) )**((k+1)/(k-1)))**0.5) / ((P_exit_3/P_max_3)**(1.0/k) * (((2.0*k) / (k-1)) * (1.0 - (P_exit_3/P_max_3) ** ((k-1)/k)) )**0.5))
    D_exit_3 = ((A_exit_3/math.pi)**0.5)*2.0

    print("Motor 3 throat diameter: {0} mm".format(D_throat_3*1000.0))
    print("Motor 3 exit diameter: {0} mm".format(D_exit_3*1000.0))
    """

    P_c_1 = P_max_1
    #P_c_2 = P_max_2
    #P_c_3 = P_max_3

    h_vehicle = 0.0
    v_vehicle = 0.0

    m_vehicle = m_total

    rho_1 = 1.22
    P_1 = P_exit_1
    T_1 = 273.0 + 20.0

    #rho_2 = 0.98
    #P_2 = P_exit_2
    #T_2 = 273.0 + 20.0

    #rho_3 = 0.66
    #P_3 = P_exit_3
    #T_3 = 273.0 + 20.0


    while rod_stage_1 > 5.0/1000.0:

        P_c_1 = (Decimal(A_burn_intermediate_1)/Decimal(A_throat_1)**(Decimal(1.0)/(Decimal(1.0) - Decimal(n_motor))))
        P_c_intermediate_1.append(float(P_c_1)/1000000.0)

        r_motor_1 = (burn_rate(P_c_1, a_motor, n_motor))
        r_motor_intermediate_1.append(r_motor_1*1000.0)

        m_flow_intermediate_1.append(A_burn_intermediate_1*r_motor_1*rho_grain_NG_NC)
        m_flow_total = m_flow_total + A_burn_intermediate_1*r_motor_1*rho_grain_NG_NC*dt

        C_f_1 = ((2.0*k**2.0)/(k-1.0)*(2.0/(k+1.0))**((k+1)/(k-1))*(1.0-(P_exit_1/float(P_c_1))**((k-1.0)/k)))**0.5 + ((P_exit_1 - P_1)/float(P_c_1)) * (A_exit_1/A_throat_1)
        F_thrust_intermediate_1.append(C_f_1*A_throat_1*float(P_c_1))

        A_burn_intermediate_1 = (2.0*math.pi*(R_rod - r_motor_1*n_counter/n_steps) + 2.0*math.pi*(R_tube + r_motor_1*n_counter/n_steps)) * L_grain_1
        rod_stage_1 = R_rod - r_motor_1*n_counter/n_steps
        erosive_factor_1.append((math.pi*(R_tube + r_motor_1*n_counter/n_steps)**2.0 - math.pi*(R_rod - r_motor_1*n_counter/n_steps)**2.0)/A_throat_1)

        rho_1,P_1,T_1 = compute_density(h_vehicle)
        a_vehicle = (C_f_1*A_throat_1*float(P_c_1) - get_gravity(h_vehicle)*m_vehicle - (D_capsule/2)**2*math.pi*Cd_z*rho_1*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - (A_burn_intermediate_1*r_motor_1*rho_grain_NG_NC)*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        rho_exhaust_1 = float(P_c_1)/(R_gas_specific*T_flame)
        M_exhaust_1.append(((m_flow_1/rho_exhaust_1)/(math.pi*((D_max/2.0)**2.0)))/(k*R_gas_specific*T_flame)**0.5)

        Isp_stage_1.append((C_f_1*A_throat_1*float(P_c_1)) / (get_gravity(h_vehicle) * (A_burn_intermediate_1*r_motor_1*rho_grain_NG_NC)))


        n_counter = n_counter + 1.0
        timestamp_1.append(t_burn)
        t_burn = t_burn + dt



    plt.figure()
    plt.plot(timestamp_1, F_thrust_intermediate_1)
    plt.title("Stage 1 Thrust")
    plt.xlabel("Time (s)")
    plt.ylabel("Thrust (N)")

    plt.figure()
    plt.plot(timestamp_1, P_c_intermediate_1, label='Pressure')
    plt.plot(timestamp_1, r_motor_intermediate_1, label='Burn rate')
    plt.plot(timestamp_1, M_exhaust_1, label='Input mach number')
    plt.plot(timestamp_1, erosive_factor_1, label='Erosion burning factor')


    plt.legend()

    plt.title("Stage 1 Chamber Pressure & burn rate")
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure (MPa)/Burn Rate (mm/s)")

    plt.figure()
    plt.plot(timestamp_1, m_flow_intermediate_1)
    plt.title("Stage 1 Massflow")
    plt.xlabel("Time (s)")
    plt.ylabel("Massflow (kg/s)\n")

    plt.figure()
    plt.plot(timestamp_1, Isp_stage_1)
    plt.title("Stage 1 Specific Impulse")
    plt.xlabel("Time (s)")
    plt.ylabel("Isp (s)\n")
    
    
    #m_vehicle = m_vehicle - m_individual_stage_structural_1

    print("Rho 1: {0} kg/m3".format(rho_1))

    """

    t_burn = 0.0
    n_counter = 0.0
    altitude = 0.0

    while rod_stage_2 > 5.0/1000.0:

        P_c_2 = (Decimal(A_burn_intermediate_2)/Decimal(A_throat_2)**(Decimal(1.0)/(Decimal(1.0) - Decimal(n_motor))))
        P_c_intermediate_2.append(float(P_c_2)/1000000.0)

        r_motor_2 = (burn_rate(P_c_2, a_motor, n_motor))
        r_motor_intermediate_2.append(r_motor_2*1000.0)

        m_flow_intermediate_2.append(A_burn_intermediate_2*r_motor_2*rho_grain_NG_NC)
        m_flow_total = m_flow_total + A_burn_intermediate_2*r_motor_2*rho_grain_NG_NC*dt

        C_f_2 = ((2.0*k**2.0)/(k-1.0)*(2.0/(k+1.0))**((k+1)/(k-1))*(1.0-(P_exit_2/float(P_c_2))**((k-1.0)/k)))**0.5 + ((P_exit_2 - P_2)/float(P_c_2)) * (A_exit_2/A_throat_2)
        F_thrust_intermediate_2.append(C_f_2*A_throat_2*float(P_c_2))

        A_burn_intermediate_2 = (2.0*math.pi*(R_rod - r_motor_2*n_counter/n_steps) + 2.0*math.pi*(R_tube + r_motor_2*n_counter/n_steps)) * L_grain_2
        rod_stage_2 = R_rod - r_motor_2*n_counter/n_steps

        erosive_factor_2.append((math.pi*(R_tube + r_motor_2*n_counter/n_steps)**2.0 - math.pi*(R_rod - r_motor_2*n_counter/n_steps)**2.0)/A_throat_2)

        rho_2,P_2,T_2 = compute_density(h_vehicle)
        a_vehicle = (C_f_2*A_throat_2*float(P_c_2) - get_gravity(h_vehicle)*m_vehicle - (D_capsule/2.0)**2*math.pi*Cd_z*rho_2*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - (A_burn_intermediate_2*r_motor_2*rho_grain_NG_NC)*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        rho_exhaust_2 = float(P_c_2)/(R_gas_specific*T_flame)
        M_exhaust_2.append(((m_flow_2/rho_exhaust_2)/(math.pi*((D_max/2.0)**2.0)))/(k*R_gas_specific*T_flame)**0.5)

        n_counter = n_counter + 1.0
        timestamp_2.append(t_burn)
        t_burn = t_burn + dt

   
    plt.figure()
    plt.plot(timestamp_2, F_thrust_intermediate_2)
    plt.title("Stage 2 Thrust")
    plt.xlabel("Time (s)")
    plt.ylabel("Thrust (N)")

    plt.figure()
    plt.plot(timestamp_2, P_c_intermediate_2, label='Pressure')
    plt.plot(timestamp_2, r_motor_intermediate_2, label='Burn rate')
    plt.plot(timestamp_2, M_exhaust_2, label='Input mach number')
    plt.plot(timestamp_2, erosive_factor_2, label='Erosion burning factor')


    plt.legend()

    plt.title("Stage 2 Chamber Pressure & burn rate")
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure (MPa)/Burn Rate (mm/s)")

    plt.figure()
    plt.plot(timestamp_2, m_flow_intermediate_2)
    plt.title("Stage 2 Massflow")
    plt.xlabel("Time (s)")
    plt.ylabel("Massflow (kg/s)\n")


    m_vehicle = m_vehicle - m_individual_stage_structural_2

    t_burn = 0.0
    n_counter = 0.0
    altitude = 0.0


    print("Rho 2: {0} kg/m3".format(rho_2))

    while rod_stage_3 > 5.0/1000.0:

        r_motor_3 = (burn_rate(P_c_3, a_motor, n_motor))
        r_motor_intermediate_3.append(r_motor_3*1000.0)

        m_flow_intermediate_3.append(A_burn_intermediate_3*r_motor_3*rho_grain_NG_NC)
        m_flow_total = m_flow_total + A_burn_intermediate_3*r_motor_3*rho_grain_NG_NC*dt

        C_f_3 = ((2*k**2.0)/(k-1.0)*(2.0/(k+1.0))**((k+1)/(k-1))*(1.0-(P_exit_3/float(P_c_3))**((k-1.0)/k)))**0.5 + ((P_exit_3 - P_3)/float(P_c_3)) * (A_exit_3/A_throat_3)
        F_thrust_intermediate_3.append(C_f_3*A_throat_3*float(P_c_3))

        A_burn_intermediate_3 = (2.0*math.pi*(R_rod - r_motor_3*n_counter/n_steps) + 2.0*math.pi*(R_tube + r_motor_3*n_counter/n_steps)) * L_grain_3
        rod_stage_3 = R_rod - r_motor_3*n_counter/n_steps

        erosive_factor_3.append((math.pi*(R_tube + r_motor_3*n_counter/n_steps)**2.0 - math.pi*(R_rod - r_motor_3*n_counter/n_steps)**2.0)/A_throat_3)

        P_c_3 = (Decimal(A_burn_intermediate_3)/Decimal(A_throat_3)**(Decimal(1.0)/(Decimal(1.0) - Decimal(n_motor))))
        P_c_intermediate_3.append(float(P_c_3)/1000000.0)

        rho_3,P_3,T_3 = compute_density(h_vehicle)
        a_vehicle = (C_f_3*A_throat_3*float(P_c_3) - get_gravity(h_vehicle)*m_vehicle - (D_capsule/2)**2*math.pi*Cd_z*rho_3*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - (A_burn_intermediate_3*r_motor_3*rho_grain_NG_NC)*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        rho_exhaust_3 = float(P_c_3)/(R_gas_specific*T_flame)
        M_exhaust_3.append(((m_flow_3/rho_exhaust_3)/(math.pi*((D_max/2.0)**2.0)))/(k*R_gas_specific*T_flame)**0.5)

        n_counter = n_counter + 1.0
        timestamp_3.append(t_burn)
        t_burn = t_burn + dt


    plt.figure()
    plt.plot(timestamp_3, F_thrust_intermediate_3)
    plt.title("Stage 3 Thrust")
    plt.xlabel("Time (s)")
    plt.ylabel("Thrust (N)")

    plt.figure()
    plt.plot(timestamp_3, P_c_intermediate_3, label='Pressure')
    plt.plot(timestamp_3, r_motor_intermediate_3, label='Burn rate')
    plt.plot(timestamp_3, M_exhaust_3, label='Input mach number')
    plt.plot(timestamp_3, erosive_factor_3, label='Erosion burning factor')


    plt.legend()

    plt.title("Stage 3 Chamber Pressure & burn rate")
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure (MPa)/Burn Rate (mm/s)")

    plt.figure()
    plt.plot(timestamp_3, m_flow_intermediate_3)
    plt.title("Stage 3 Massflow")
    plt.xlabel("Time (s)")
    plt.ylabel("Massflow (kg/s)\n")

    plt.show()
    """


    print("The burn rate of the motor at {0} MPa is {1} mm/s".format(P_max_1/(1.0e6), r_motor_1*1000.0))
    print("Number of stages: {0}".format(n_stages))
    
    print("Burn time 1 stage: {0} seconds\n".format(t_burn))
    #print("Burn time 2 stage: {0} seconds\n".format(t_stage_2))
    #print("Burn time 3 stage: {0} seconds\n".format(t_stage_3))

    print("Massflow stage 1: {0} kg/s".format(m_flow_1))
    #print("Massflow stage 2: {0} kg/s".format(m_flow_2))
    #print("Massflow stage 3: {0} kg/s\n".format(m_flow_3))

    print("Burn surface stage 1: {0} m2".format(A_burn_1))
    #print("Burn surface stage 2: {0} m2".format(A_burn_2))
    #print("Burn surface stage 3: {0} m2\n".format(A_burn_3))

    print("Grain length stage 1: {0} mm".format(L_grain_1*1000.0))
    #print("Grain length stage 2: {0} mm".format(L_grain_2*1000.0))
    #print("Grain length stage 3: {0} mm\n".format(L_grain_3*1000.0))

    print("Starting diameter rod stage 1-2: {0} mm".format(R_rod*2.0*1000.0))
    print("Starting diameter tube stage 1-2: {0} mm".format(R_tube*2.0*1000.0))
    print("Final diameter tube stage 1-2: {0} mm".format(motor_OD_final*1000.0))
    
    print("Starting rod-tube gap stage 1-2: {0} mm".format(W_tube*1000.0))

    print("Specific gas constant: {0} J/k.kg".format(R_gas_specific))
    print("Cp heat capacity: {0} J/k.kg".format(Cp))
    print("Molar mass: {0} kg/mol".format(M_gas))

    print("Nozzle 1 length: {0} mm".format(435.0))
    #print("Nozzle 2 length: {0} mm".format(421.0))
    #print("Nozzle 3 length: {0} mm\n".format(494.0))

    print("Total burned propellant mass: {0} kg".format(m_flow_total))
    print("Total loaded propellant mass: {0} kg\n".format(m_grain_total))
    
    plt.show()


def number_of_stages(burn_rate, R_rod, t_total, R_rod_support):
    
    burn_rate = float(burn_rate)
    t_total = float(t_total)
    
    t_stage = (R_rod-R_rod_support)/burn_rate
    return t_total/t_stage, t_stage

def burn_rate(P_motor, a_motor, n_motor):

    a_motor = float(a_motor)
    n_motor = float(n_motor)
    P_motor = float(P_motor)

    return a_motor*P_motor**n_motor

def compute_density(altitude):

    T0 = 293
    P0 = 101250.0
    Rho_0 = 1.225

    R_ideal = 8.31447
    M_air = 0.0289644
    g0 = 9.81

    L = 0.0065
    T_altitude = T0 - L*altitude
    P_altitude = P0 * (1 - (L*altitude)/T0) ** ((g0 * M_air) / (R_ideal * L))
    rho_altitude = (P_altitude * M_air) / (R_ideal * T_altitude)

    return rho_altitude, P_altitude, T_altitude

def get_gravity(altitude):

    G_grav = 6.67408 * 10**-11
    M_earth = 5.972*10**24
    R_earth = 6371.0*10**3

    return (G_grav * M_earth)/((R_earth + altitude)**2)



main()