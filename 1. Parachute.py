import math
import matplotlib.pyplot as plt
import numpy as np
from decimal import *

def main():

    #Parachute and airbrake design
    h_burnout_stage_1 = 8441.86
    h_burnout_stage_2 = 10549.0 
    h_burnout_stage_3 = 28859.0 

    v_burnout_stage_1 = 1670.0
    v_burnout_stage_2 = 839.0 
    v_burnout_stage_3 = 1592.0 

    h_apogee_payload = 117.0*10**3

    h_apogee_stage_1 = h_burnout_stage_1
    h_apogee_stage_2 = h_burnout_stage_2
    h_apogee_stage_3 = h_burnout_stage_3

    v_descent_1 = [0.0]
    v_descent_2 = [0.0]
    v_descent_3 = [0.0]
    

    m_payload = 50.0 

    m_individual_stage_structural_1 = 170.0
    m_individual_stage_structural_2 = 27.5
    m_individual_stage_structural_3 = 16.5

    Cd_parachute = 1.25
    Cd_stage = 0.03
    D_stage = 0.442

    v_touchdown = 10.0 #m/s

    dt = 1.0/1000.0
    g0 = 9.81
    a_acc = 0.0
    a_acc_down = 0.0



    while v_burnout_stage_1 > 0.0:
        
        rho_altitude_1, P_altitude_1, T_altitude_1 = compute_density(h_apogee_stage_1)

        a_acc = ((0.0 - m_individual_stage_structural_1*get_gravity(h_apogee_stage_1) -
                                                                    Cd_stage*(math.pi*(D_stage/2.0)**2.0)*0.5*rho_altitude_1*v_burnout_stage_1**2.0))/m_individual_stage_structural_1
        
        v_burnout_stage_1 = v_burnout_stage_1 + a_acc*dt
        h_apogee_stage_1 = h_apogee_stage_1 + v_burnout_stage_1*dt


    print("Stage 1 apogee altitude: {0} m".format(h_apogee_stage_1))

    """
    while v_burnout_stage_2 > 0.0:
        
        rho_altitude_2, P_altitude_2, T_altitude_2 = compute_density(h_apogee_stage_2)

        a_acc = ((0.0 - m_individual_stage_structural_2*get_gravity(h_apogee_stage_2) -
                                                                    Cd_stage*(math.pi*(D_stage/2.0)**2.0)*0.5*rho_altitude_2*v_burnout_stage_2**2.0))/m_individual_stage_structural_2
        
        v_burnout_stage_2 = v_burnout_stage_2 + a_acc*dt
        h_apogee_stage_2 = h_apogee_stage_2 + v_burnout_stage_2*dt


    print("Stage 2 apogee altitude: {0} m".format(h_apogee_stage_2))


    while v_burnout_stage_3 > 0.0:
        
        rho_altitude_3, P_altitude_3, T_altitude_3 = compute_density(h_apogee_stage_3)

        a_acc = ((0.0 - m_individual_stage_structural_3*get_gravity(h_apogee_stage_3) -
                                                                    Cd_stage*(math.pi*(D_stage/2.0)**2.0)*0.5*rho_altitude_3*v_burnout_stage_3**2.0))/m_individual_stage_structural_3
        
        v_burnout_stage_3 = v_burnout_stage_3 + a_acc*dt
        h_apogee_stage_3 = h_apogee_stage_3 + v_burnout_stage_3*dt


    print("Stage 3 apogee altitude: {0} m\n".format(h_apogee_stage_3))
    """

    #Compute the descent velocity for each parachute
    D_parachute_1 = 3000.0/1000.0
    D_parachute_2 = 3000.0/1000.0
    D_parachute_3 = 3000.0/1000.0

    number_parachutes_1 = 2.0

    t_descent_1 = [0.0]
    t_descent_2 = [0.0]
    t_descent_3 = [0.0]

    counter = 0

    while h_apogee_stage_1 > 0.0:
        
        rho_altitude_1, P_altitude_1, T_altitude_1 = compute_density(h_apogee_stage_1)
        A_parachute_1 = math.pi*(D_parachute_1/2.0)**2.0

        a_acc_down = (get_gravity(h_apogee_stage_1)*m_individual_stage_structural_1 - number_parachutes_1*(A_parachute_1 * Cd_parachute*0.5*rho_altitude_1*v_burnout_stage_1**2.0))/m_individual_stage_structural_1
        v_burnout_stage_1 = v_burnout_stage_1 + a_acc_down*dt
        h_apogee_stage_1 = h_apogee_stage_1 - v_burnout_stage_1*dt

        v_descent_1.append(v_burnout_stage_1)
        t_descent_1.append(t_descent_1[counter] + dt)

        counter = counter + 1

    plt.figure()
    plt.plot(t_descent_1, v_descent_1)
    plt.title("Stage 1 descent velocity")
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")

    """

    counter = 0

    while h_apogee_stage_2 > 0.0:
        
        rho_altitude_2, P_altitude_2, T_altitude_2 = compute_density(h_apogee_stage_2)
        A_parachute_2 = math.pi*(D_parachute_2/2.0)**2.0

        a_acc_down = (get_gravity(h_apogee_stage_2)*m_individual_stage_structural_2 - A_parachute_2 * Cd_parachute*0.5*rho_altitude_2*v_burnout_stage_2**2.0)/m_individual_stage_structural_2
        v_burnout_stage_2 = v_burnout_stage_2 + a_acc_down*dt
        h_apogee_stage_2 = h_apogee_stage_2 - v_burnout_stage_2*dt

        v_descent_2.append(v_burnout_stage_2)
        t_descent_2.append(t_descent_2[counter] + dt)

        counter = counter + 1
    
    plt.figure()
    plt.plot(t_descent_2, v_descent_2)
    plt.title("Stage 2 descent velocity")
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")




    counter = 0

    while h_apogee_stage_3 > 0.0:
        
        rho_altitude_3, P_altitude_3, T_altitude_3 = compute_density(h_apogee_stage_3)
        A_parachute_3 = math.pi*(D_parachute_3/2.0)**2.0

        a_acc_down = (get_gravity(h_apogee_stage_3)*m_individual_stage_structural_3 - A_parachute_3 * Cd_parachute*0.5*rho_altitude_3*v_burnout_stage_3**2.0)/m_individual_stage_structural_3
        v_burnout_stage_3 = v_burnout_stage_3 + a_acc_down*dt
        h_apogee_stage_3 = h_apogee_stage_3 - v_burnout_stage_3*dt

        v_descent_3.append(v_burnout_stage_3)
        t_descent_3.append(t_descent_3[counter] + dt)

        counter = counter + 1

    plt.figure()
    plt.plot(t_descent_3, v_descent_3)
    plt.title("Stage 3 descent velocity")
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")

    """

    plt.show()


def compute_density(altitude):

    T0 = 293.0
    P0 = 101250.0
    Rho_0 = 1.225

    R_ideal = 8.31447
    M_air = 0.0289644
    g0 = 9.81

    L = 0.0065
    T_altitude = T0 - L*altitude
    
    if ((1.0 - (L*altitude)/T0)) > 0.0000000001:

        P_altitude = P0 * (1.0 - (L*altitude)/T0) ** ((g0 * M_air) / (R_ideal * L))
        rho_altitude = (P_altitude * M_air) / (R_ideal * T_altitude)

    else:
        P_altitude = 0.0
        rho_altitude = 0.0

    return rho_altitude, P_altitude, T_altitude


def get_gravity(altitude):

    G_grav = 6.67408 * 10**-11
    M_earth = 5.972*10**24
    R_earth = 6371.0*10**3

    return (G_grav * M_earth)/((R_earth + altitude)**2)



main()