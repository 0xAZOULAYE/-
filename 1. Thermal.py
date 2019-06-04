import math
import matplotlib.pyplot as plt
import numpy as np
from decimal import *


def main():

    #This script is meant to dimension the thickness of the motor heat shield

    T_comb = 3070.0 #Motor combustion temperature in K
    D_core_motor = 0.4
    D_outer_airframe = 1.0 #m

    T_max_airframe = 293.0 #K maximum allowable temperature of the airframe (Al-7075)

    Mat_airframe = "7075-T6"
    Mat_casing = "304"
    Mat_heat_shield = "Fiberglass"
    Mat_grain = "NC + NG + DEP + 2NDPA"

    k_airframe_material = 0.0 #Heat transfer coefficient of the airframe material
    k_motor_casing_material = 0.0 #Heat transfer coefficient of the motor casing material
    k_motor_heat_shield_material = 0.0 #Heat transfer coefficient of the motor heat shield material
    k_motor_grain_material = 0.0 #Heat transfer coefficient of the motor grain

    #Define the motor operational characteristics
    t_burn = 12.0 #Seconds

    M_gas = 0.0
    R_ideal = 8.31447 # J/mol*K
    g = 9.81 #m/s

    #Compute the motor heating
    T_comb = 2760.0
    H_exp = 4.36e6 #J/kg

    #Compute the aerodynamic heating---------------------------------------------------------------------------------------------------------------
    #Page 930 of fundamentals of aerodynamics

    v_body = [0.0]
    T_body = [0.0]
    timestamp = [0.0]

    alpha_max = 10.0/(180/math.pi)
    Cd_body_side = 1.143
    Cl_body_side = 0.382
    Cd_z = 0.03
    r_body = 442.0/2000.0
    r_core = 432.0/2000.0

    F_thrust_1 = 52000.0 #N
    #F_thrust_2 = 21000.0 #N
    #F_thrust_3 = 12600.0 #N

    massflow_motor_1 = 30.3 #kg/s 
    #massflow_motor_2 = 7.24
    #massflow_motor_3 = 3.14

    m_individual_stage_structural_1 = 170.0
    #m_individual_stage_structural_2 = 27.5
    #m_individual_stage_structural_3 = 16.5

    m_total = 561.0
    m_vehicle = m_total

    t_flight = 0.0
        
    t_burn_1 = 12.0
    #t_burn_2 = 14.0
    #t_burn_3 = 16.0
        
    dt = 1.0/1000.0

    a_vehicle = 0.0
    v_vehicle = 0.0
    h_vehicle = 0.0

    T_stagnation_cone = [0.0]
    gamma = 1.2
    R_air = 287.0

    while (t_flight < t_burn_1):

        rho_1,P_1,T_1 = compute_density(h_vehicle)
        a_vehicle = (F_thrust_1 - 9.81*m_vehicle - (0.442/2.0)**2.0*math.pi*Cd_z*rho_1*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - massflow_motor_1*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        v_body.append(v_vehicle)
        timestamp.append(t_flight)

        v_sound = (gamma*R_air*T_1)**0.5
        M_cone = v_vehicle/v_sound
        T_stagnation_cone.append((T_1 * (1.0 + (gamma-1.0)/2.0*M_cone**2.0)) - 273.0)

        t_flight = t_flight + dt

    m_vehicle = m_vehicle - m_individual_stage_structural_1
    """

    while (t_flight < t_burn_2):

        rho_2,P_2,T_2 = compute_density(h_vehicle)
        a_vehicle = (F_thrust_2 - 9.81*m_vehicle - (0.3/2.0)**2.0*math.pi*Cd_z*rho_2*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - massflow_motor_2*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        v_body.append(v_vehicle)
        timestamp.append(t_flight)

        v_sound = (gamma*R_air*T_2)**0.5
        M_cone = v_vehicle/v_sound
        T_stagnation_cone.append((T_2 * (1.0 + (gamma-1.0)/2.0*M_cone**2.0)) - 273.0)

        t_flight = t_flight + dt


    m_vehicle = m_vehicle - m_individual_stage_structural_1 - m_individual_stage_structural_2
    

    while (t_flight < t_burn_3):

        rho_3,P_3,T_3 = compute_density(h_vehicle)
        a_vehicle = (F_thrust_3 - 9.81*m_vehicle - (0.3/2.0)**2.0*math.pi*Cd_z*rho_3*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - massflow_motor_3*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        v_body.append(v_vehicle)
        timestamp.append(t_flight)

        v_sound = (gamma*R_air*T_3)**0.5
        M_cone = v_vehicle/v_sound
        T_stagnation_cone.append((T_3 * (1.0 + (gamma-1.0)/2.0*M_cone**2.0)) - 273.0)

        t_flight = t_flight + dt
    """


    plt.figure()
    plt.plot(timestamp, v_body)
    plt.title("Rocket velocity")
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")

    plt.figure()
    plt.plot(timestamp, T_stagnation_cone)
    plt.title("Rocket cone temperature")
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (C)")

    plt.show()


    #Compute the body temperature
    T_body_0 = 293.0
    T_body_1 = 293.0 + 50.0

    mat_mech_name = "7075-T6"
    mat_mech_E = 69.0e9
    mat_mech_sigma_ys = 430.0e6
    mat_mech_sigma_shear = mat_mech_sigma_ys*1.0
    mat_mech_sigma_fat = 160.0e6
    mat_mech_max_T = 293.0 #K
    mat_mech_CT = 960.0 #Thermal capacity
    mat_mech_rho = 2810.0

    mat_nozzle_name = "graphite"
    mat_nozzle_E = 8.0e9
    mat_nozzle_sigma_ys = 6.9e6
    mat_nozzle_sigma_shear = mat_nozzle_sigma_ys*0.6
    mat_nozzle_sigma_fat = 0.0e6
    mat_nozzle_max_T = 0.0 #K
    mat_nozzle_CT = 710.0 #Thermal capacity
    mat_nozzle_k = 25.0 #w/m.k
    
    mat_insulation_name = "fiberglass"
    mat_insulation_E = 0.0
    mat_insulation_sigma_ys = 1080.0e6
    mat_insulation_sigma_shear = mat_insulation_sigma_ys*0.6
    mat_insulation_sigma_fat = 0.0
    mat_insulation_max_T = 1000.0 #K
    mat_insulation_CT = 700.0 #Thermal capacity
    mat_insulation_rho = 144.0
    mat_insulation_k = 0.036 #W/m.k

    L_stage_1 = 5.0 #2692.0/1000.0
    #L_stage_2 = 5.0*340.0/1000.0 #1692.0/1000.0
    #L_stage_3 = 4.0*340.0/1000.0 #1308.0/1000.0
    
    t_wall = 1.75/1000.0

    V_stage_1 = ((0.432/2.0 + t_wall)**2.0*math.pi - (0.432/2.0)**2.0*math.pi)*L_stage_1
    #V_stage_2 = ((0.35/2.0 + t_wall)**2.0*math.pi - (0.35/2.0)**2.0*math.pi)*L_stage_2
    #V_stage_3 = ((0.35/2.0 + t_wall)**2.0*math.pi - (0.35/2.0)**2.0*math.pi)*L_stage_3

    Q_stage_1 = mat_mech_CT*V_stage_1*mat_mech_rho*(T_body_1-T_body_0)
    #Q_stage_2 = mat_mech_CT*V_stage_2*mat_mech_rho*(T_body_1-T_body_0)
    #Q_stage_3 = mat_mech_CT*V_stage_3*mat_mech_rho*(T_body_1-T_body_0)

    q_stage_1 = Q_stage_1/t_burn_1
    #q_stage_2 = Q_stage_2/t_burn_2
    #q_stage_3 = Q_stage_3/t_burn_3

    R_insulation_1 = 0.0
    #R_insulation_2 = 0.0
    #R_insulation_3 = 0.0

    R_insulation_1 = math.exp((2.0*math.pi*L_stage_1*mat_insulation_k*(T_comb - 293.0))/q_stage_1)*r_core
    #R_insulation_2 = math.exp((2.0*math.pi*L_stage_2*mat_insulation_k*(T_comb - 293.0))/q_stage_2)*r_core
    #R_insulation_3 = math.exp((2.0*math.pi*L_stage_3*mat_insulation_k*(T_comb - 293.0))/q_stage_3)*r_core

    t_insulation_1 =  R_insulation_1 - r_core
    #t_insulation_2 =  R_insulation_2 - r_core
    #t_insulation_3 =  R_insulation_3 - r_core

    print("Stage 1 insulation thickness: {0} mm".format(t_insulation_1*1000.0))
    #print("Stage 2 insulation thickness: {0} mm".format(t_insulation_2*1000.0))
    #print("Stage 3 insulation thickness: {0} mm".format(t_insulation_3*1000.0))

    m_insulation_1 = ((math.pi*R_insulation_1**2.0 - math.pi*r_core**2.0)*L_stage_1 + math.pi*R_insulation_1**2.0*t_insulation_1*2.0)*mat_insulation_rho
    #m_insulation_2 = ((math.pi*R_insulation_2**2.0 - math.pi*r_core**2.0)*L_stage_2 + math.pi*R_insulation_2**2.0*t_insulation_2*2.0)*mat_insulation_rho
    #m_insulation_3 = ((math.pi*R_insulation_3**2.0 - math.pi*r_core**2.0)*L_stage_3 + math.pi*R_insulation_3**2.0*t_insulation_3*2.0)*mat_insulation_rho

    print("Insulation mass stage 1: {0} kg".format(m_insulation_1))
    #print("Insulation mass stage 2: {0} kg".format(m_insulation_2))
    #print("Insulation mass stage 3: {0} kg".format(m_insulation_3))


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

main()





















