import math
import matplotlib.pyplot as plt
import numpy as np
from decimal import *

def main():

    alpha_max = 10.0/(180/math.pi)
    Cd_body_side = 1.143
    Cl_body_side = 0.382
    Cd_z = 0.03
    r_body = 0.442/2.0

    dynamic_pressure = [0.0]

    F_thrust_1 = 52000.0 #N
    #F_thrust_2 = 21000.0 #N
    #F_thrust_3 = 12600.0 #N

    massflow_motor_1 = 30.3 #kg/s 
    #massflow_motor_2 = 7.24
    #massflow_motor_3 = 3.14

    v_wind = 10.0

    #Calculate center of mass
    m_pay = 50.0
    d_pay = 5421.0/1000.0

    m_prop_1 = 370.0
    d_prop_1 = 2491.0/1000.0

    #m_prop_2 = (350.0/(340.0*3 + 340.0*5 + 340.0*8.0)) * (340.0*5.0)
    #d_prop_2 = 4.726

    #m_prop_3 = (350.0/(340.0*3 + 340.0*5 + 340.0*8.0)) * (340.0*4.0)
    #d_prop_3 = 6.997

    m_individual_stage_structural_1 = 561.0 - m_prop_1 - m_pay
    #m_individual_stage_structural_2 = 27.5
    #m_individual_stage_structural_3 = 16.5
    

    m_airframe = m_individual_stage_structural_1 #+ m_individual_stage_structural_2 + m_individual_stage_structural_3
    d_airframe = 5421.0/2000.0

    m_total = m_individual_stage_structural_1 + 350.0 + m_pay #+ m_individual_stage_structural_2 + m_individual_stage_structural_3
    d_mass = ((m_prop_1*d_prop_1) + m_pay*d_pay + m_airframe*d_airframe) / m_total

    #Calculate center of pressure

    A_n = 1.0*0.442*0.5
    d_n = 5898.0/1000.0

    A_b = 0.442*5.421
    d_b = 5421.0/2000.0

    A_f = (0.820-0.5) * 0.5 * 2.0
    d_f = (0.820-0.5)/2.0

    #A_f2 = 0.34*0.16*2.0
    #d_f2 = 4.062

    #A_f3 = 0.34*0.16*2.0
    #d_f3 = 6.472

    A_p = A_n + A_b + A_f
    d_p = ((A_n*d_n) + (A_b*d_b) + (A_f*d_f)) / A_p

    Cd_fin = 1.28
    t_fin = 5.0/1000.0

    #Plot the first stage forces and angles of attack
    t_flight = 0.0
    
    t_burn_1 = 12.2
    #t_burn_2 = 12.09
    #t_burn_3 = 16.77
    
    dt = 1.0/1000.0

    h_vehicle = 0.0
    v_vehicle = 0.0
    a_vehicle = 0.0
    m_vehicle = m_total

    L_stage_1 = 5463.0/1000.0
    #L_stage_2 = 5072.0/1000.0
    #L_stage_3 = 2662.0/1000.0

    I_stage_1_end = 1.0/3.0 * m_total * L_stage_1**2.0
    I_stage_1_middle = 1.0/12.0 * m_total * L_stage_1**2.0
    I_stage_1_rot = 0.5*m_total*0.15**2.0

    alpha_vehicle = [0.0]
    timestamp_1 = [0.0]
    counter = 0

    angular_acc_1 = 0.0
    angular_velo_1 = 0.0

    angle_momentary_1 = 0.0
    velo_momentary_1 = 0.0
    acc_momentary_1 = 0.0

    L_downrange = [0.0]

    roughness = 12.5*10**-6.0
    thickness = 5.0*10**-3.0

    omega_rocket = 100.0 #RPS
    omega_rocket_rad = 100.0*2.0*math.pi #rad/s
    alpha_airfoil = 1.0 #max alpha = 12 degrees

    l_rotation_fin = 150.0/1000.0
    w_rotation_fin = 10.0/1000.0
    number_rotation_fin = 2.0

    A_rotation_fin = l_rotation_fin*w_rotation_fin*number_rotation_fin

    cl_flat = 2.0*math.pi*math.sin(alpha_airfoil/(180.0/math.pi))
    cd_flat = cl_flat + 0.43/((math.log(100.0/roughness))**2.56) + 0.3*thickness*math.cos(alpha_airfoil)

    rot_acc = 0.0
    rot_velo = [0.0]

    t_stage = 0.0

    print("The center of pressure is located {0} meters from the nozzle exit (Stage 1)".format(d_p))
    print("The center of mass is located {0} meters from the nozzle exit (Stage 1)".format(d_mass))

    #First stage
    while (t_stage < t_burn_1):

        rho_1,P_1,T_1 = compute_density(h_vehicle)
        a_vehicle = (F_thrust_1 - get_gravity(h_vehicle)*m_vehicle - 0.5*(0.32/2.0)**2.0*math.pi*Cd_z*rho_1*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - massflow_motor_1*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        M_g_1 = math.sin(alpha_vehicle[counter]/(180.0/math.pi)) * get_gravity(h_vehicle) * d_mass

        L_p_1 = d_mass - d_p
        F_p_wind_1 = 0.5 * A_p * Cd_body_side * rho_1 * v_wind**2.0
        M_p_wind_1 = F_p_wind_1 * L_p_1

        acc_momentary_1 = (M_p_wind_1)/I_stage_1_middle + M_g_1/I_stage_1_end
        velo_momentary_1 = acc_momentary_1*dt
        angle_momentary_1 = velo_momentary_1*dt

        F_lift_1 = 0.5*rho_1*Cl_body_side*(v_vehicle*math.sin(angle_momentary_1))**2.0*A_p
        F_drag_1 = 0.5*rho_1*Cd_body_side*(v_vehicle*math.sin(angle_momentary_1))**2.0*A_p

        F_ld_1 = -1.0 * (F_lift_1**2.0 + F_drag_1**2.0)**0.5
        M_ld_1 = (d_mass - d_p) * F_ld_1

        angular_acc_1 = (M_ld_1)/I_stage_1_middle
        angular_velo_1 = velo_momentary_1 + angular_acc_1*dt

        alpha_vehicle.append(alpha_vehicle[counter - 1] + angular_velo_1*dt*(180.0/math.pi))
        L_downrange.append(math.cos((alpha_vehicle[counter - 1] + angular_velo_1*dt*(180.0/math.pi))/(180.0/math.pi)) * v_vehicle * 1.0/1000.0)

        m_prop_1 = (m_vehicle/(270.0*2 + 270.0*3 + 270.0*11.0)) * (270.0*11.0)
        d_mass = (((m_prop_1 - t_flight*massflow_motor_1)*d_prop_1) + m_pay*d_pay + m_airframe*d_airframe) / m_vehicle

        I_stage_1_end = 1.0/3.0 * m_vehicle * L_stage_1**2.0
        I_stage_1_middle = 1.0/12.0 * m_vehicle * L_stage_1**2.0
        I_stage_1_rot = 0.5*m_vehicle*0.15**2.0

        #Determine the rotation
        fin_force = 0.5*cl_flat*rho_1*A_rotation_fin*v_vehicle**2.0
        fin_torque = (w_rotation_fin + r_body)*fin_force

        rot_acc = fin_torque/I_stage_1_rot
        rot_velo.append(rot_velo[counter - 1] + ((rot_acc*dt)/(2.0*math.pi))*60.0)

        dynamic_pressure.append(0.5*rho_1*v_vehicle**2.0)

        counter = counter + 1
        t_flight = t_flight + dt
        t_stage = t_stage + dt
        timestamp_1.append(t_flight)


    print("\nStage 1 burnout altitude: {0} m".format(h_vehicle))
    print("Stage 1 burnout velocity: {0} m/s".format(v_vehicle))

    plt.figure()
    plt.plot(timestamp_1, rot_velo)
    plt.title("Rocket spin")
    plt.xlabel("Time (s)")
    plt.ylabel("Rotation (RPM)")

    """

    #Calculate center of mass
    m_pay = 50.0
    d_pay = 4.566

    m_prop_2 = 350.0/(340.0*3 + 340.0*5 + 340.0*8.0) * (340.0*5.0)
    d_prop_2 = 1.302

    m_prop_3 = 350.0/(340.0*3 + 340.0*5 + 340.0*8.0) * (340.0*4.0)
    d_prop_3 = 3.537

    m_airframe = m_individual_stage_structural_2 + m_individual_stage_structural_3

    d_airframe = 6.115/2.0

    m_total = 500.0
    d_mass = ((m_prop_1*d_prop_1) + (m_prop_2*d_prop_2) + (m_prop_3*d_prop_3) + m_pay*d_pay + m_airframe*d_airframe) / m_total

    m_vehicle = m_vehicle - m_individual_stage_structural_1

    A_n = 0.5*0.35*0.8
    d_n = 5.306

    A_b = 0.35*2.566
    d_b = 2.566/2.0

    A_f = (0.34*0.15) * 2.0
    d_f = 0.602

    A_f3 = 0.34*0.16
    d_f3 = 3.012

    A_p = A_n + A_b + A_f + A_f3
    d_p = ((A_n*d_n) + (A_b*d_b) + (A_f*d_f) + (A_f3*d_f3)) / A_p

    #Plot the first stage forces and angles of attack

    L_stage_1 = 8527.0/1000.0
    L_stage_2 = 5072.0/1000.0
    L_stage_3 = 2662.0/1000.0

    I_stage_2_end = 1.0/3.0 * m_total * L_stage_2**2.0
    I_stage_2_middle = 1.0/12.0 * m_total * L_stage_2**2.0

    timestamp_2 = [0.0]

    angular_acc_2 = angular_acc_1
    angular_velo_2 = angular_velo_1

    angle_momentary_2 = angle_momentary_1
    velo_momentary_2 = velo_momentary_1
    acc_momentary_2 = acc_momentary_1

    print("\nThe center of pressure is located {0} meters from the nozzle exit (Stage 2)".format(d_p))
    print("The center of mass is located {0} meters from the nozzle exit (Stage 2)".format(d_mass))

    t_stage = 0.0

    #Second stage
    while (t_stage < t_burn_2):

        rho_2,P_2,T_2 = compute_density(h_vehicle)
        a_vehicle = (F_thrust_2 - get_gravity(h_vehicle)*m_vehicle - 0.5*(0.32/2.0)**2.0*math.pi*Cd_z*rho_2*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - massflow_motor_2*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        M_g_2 = math.sin(alpha_vehicle[counter]/(180.0/math.pi)) * get_gravity(h_vehicle) * d_mass

        L_p_2 = d_mass - d_p
        F_p_wind_2 = 0.5 * A_p * Cd_body_side * rho_2 * v_wind**2.0
        M_p_wind_2 = F_p_wind_2 * L_p_2

        acc_momentary_2 = (M_p_wind_2)/I_stage_2_middle + M_g_2/I_stage_2_end
        velo_momentary_2 = acc_momentary_2*dt
        angle_momentary_2 = velo_momentary_2*dt

        F_lift_2 = 0.5*rho_2*Cl_body_side*(v_vehicle*math.sin(angle_momentary_2))**2.0*A_p
        F_drag_2 = 0.5*rho_2*Cd_body_side*(v_vehicle*math.sin(angle_momentary_2))**2.0*A_p

        F_ld_2 = -1.0 * (F_lift_2**2.0 + F_drag_2**2.0)**0.5
        M_ld_2 = (d_mass - d_p) * F_ld_2

        angular_acc_2 = (M_ld_2)/I_stage_2_middle
        angular_velo_2 = velo_momentary_2 + angular_acc_2*dt

        alpha_vehicle.append(alpha_vehicle[counter - 1] + angular_velo_2*dt*(180.0/math.pi))
        L_downrange.append(math.cos((alpha_vehicle[counter - 1] + angular_velo_2*dt*(180.0/math.pi))/(180.0/math.pi)) * v_vehicle * 1.0/1000.0)

        m_prop_2 = (m_vehicle/(270.0*2 + 270.0*3)) * (270.0*3.0)
        d_mass = (((m_prop_2 - t_flight*massflow_motor_2)*d_prop_2) + (m_prop_3*d_prop_3) + m_pay*d_pay + (m_airframe - m_individual_stage_structural_1)*d_airframe) / m_vehicle

        I_stage_2_end = 1.0/3.0 * m_vehicle * L_stage_2**2.0
        I_stage_2_middle = 1.0/12.0 * m_vehicle * L_stage_2**2.0

        dynamic_pressure.append(0.5*rho_2*v_vehicle**2.0)

        counter = counter + 1
        t_flight = t_flight + dt
        t_stage = t_stage + dt
        timestamp_1.append(t_flight)

    print("\nStage 2 burnout altitude: {0} m".format(h_vehicle))
    print("Stage 2 burnout velocity: {0} m/s".format(v_vehicle))

    #Calculate center of mass
    m_pay = 50.0
    d_pay = 2.157

    m_prop_3 = 350.0/(340.0*3 + 340.0*5 + 340.0*8.0) * (340.0*4.0)
    d_prop_3 = 1.127

    m_airframe = m_individual_stage_structural_3
    d_airframe = 2.657/2.0

    m_total = 500.0
    d_mass = ((m_prop_1*d_prop_1) + (m_prop_2*d_prop_2) + (m_prop_3*d_prop_3) + m_pay*d_pay + m_airframe*d_airframe) / m_total

    m_vehicle = m_vehicle - m_individual_stage_structural_1 - m_individual_stage_structural_2

    A_n = 0.5*0.35*0.8
    d_n = 2.895

    A_b = 0.35*2.657
    d_b = 2.657/2.0

    A_f = (0.340*0.15) * 2.0
    d_f = 0.602

    A_p = A_n + A_b + A_f
    d_p = ((A_n*d_n) + (A_b*d_b) + (A_f*d_f)) / A_p

    #Plot the first stage forces and angles of attack

    L_stage_1 = 6825.0/1000.0
    L_stage_2 = 3276.0/1000.0
    L_stage_3 = 2024.0/1000.0

    I_stage_3_end = 1.0/3.0 * m_total * L_stage_3**2.0
    I_stage_3_middle = 1.0/12.0 * m_total * L_stage_3**2.0

    timestamp_3 = [0.0]

    angular_acc_3 = angular_acc_2
    angular_velo_3 = angular_velo_2

    angle_momentary_3 = angle_momentary_2
    velo_momentary_3 = velo_momentary_2
    acc_momentary_3 = acc_momentary_2

    print("\nThe center of pressure is located {0} meters from the nozzle exit (Stage 3)".format(d_p))
    print("The center of mass is located {0} meters from the nozzle exit (Stage 3)".format(d_mass))

    t_stage = 0.0

    #Second stage
    while (t_stage < t_burn_3):

        rho_3,P_3,T_3 = compute_density(h_vehicle)
        a_vehicle = (F_thrust_3 - get_gravity(h_vehicle)*m_vehicle - 0.5*(0.32/2.0)**2.0*math.pi*Cd_z*rho_3*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - massflow_motor_3*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        M_g_3 = math.sin(alpha_vehicle[counter]/(180.0/math.pi)) * get_gravity(h_vehicle) * d_mass

        L_p_3 = d_mass - d_p
        F_p_wind_3 = 0.5 * A_p * Cd_body_side * rho_3 * v_wind**2.0
        M_p_wind_3 = F_p_wind_3 * L_p_3

        acc_momentary_3 = (M_p_wind_3)/I_stage_3_middle + M_g_3/I_stage_3_end
        velo_momentary_3 = acc_momentary_3*dt
        angle_momentary_3 = velo_momentary_3*dt

        F_lift_3 = 0.5*rho_3*Cl_body_side*(v_vehicle*math.sin(angle_momentary_3))**2.0*A_p
        F_drag_3 = 0.5*rho_3*Cd_body_side*(v_vehicle*math.sin(angle_momentary_3))**2.0*A_p

        F_ld_3 = -1.0 * (F_lift_3**2.0 + F_drag_3**2.0)**0.5
        M_ld_3 = (d_mass - d_p) * F_ld_3

        angular_acc_3 = (M_ld_3)/I_stage_3_middle
        angular_velo_3 = velo_momentary_3 + angular_acc_3*dt

        alpha_vehicle.append(alpha_vehicle[counter - 1] + angular_velo_3*dt*(180.0/math.pi))
        L_downrange.append(math.cos((alpha_vehicle[counter - 1] + angular_velo_3*dt*(180.0/math.pi))/(180.0/math.pi)) * v_vehicle * 1.0/1000.0)

        m_prop_3 = (m_vehicle/(270.0*2)) * (270.0*2.0)
        d_mass = (((m_prop_3 - t_flight*massflow_motor_3)*d_prop_3) + m_pay*d_pay + (m_airframe - m_individual_stage_structural_1 - m_individual_stage_structural_2)*d_airframe) / m_vehicle

        I_stage_3_end = 1.0/3.0 * m_vehicle * L_stage_3**2.0
        I_stage_3_middle = 1.0/12.0 * m_vehicle * L_stage_3**2.0

        dynamic_pressure.append(0.5*rho_3*v_vehicle**2.0)

        counter = counter + 1
        t_flight = t_flight + dt
        t_stage = t_stage + dt
        timestamp_1.append(t_flight)

    print("\nStage 3 burnout altitude: {0} m".format(h_vehicle))
    print("Stage 3 burnout velocity: {0} m/s".format(v_vehicle))

    """

    plt.figure()
    plt.plot(timestamp_1, alpha_vehicle)
    plt.title("Rocket angle")
    plt.xlabel("Time (s)")
    plt.ylabel("Angle (degrees)")

    plt.figure()
    plt.plot(timestamp_1, dynamic_pressure)
    plt.title("Rocket dynamic pressure")
    plt.xlabel("Time (s)")
    plt.ylabel("Q (Pa)")

    plt.figure()
    plt.plot(timestamp_1, L_downrange)
    plt.title("Rocket range")
    plt.xlabel("Time (s)")
    plt.ylabel("Range (km)")

    plt.show()



    #Define the arm length of the yo-yo despin

    omega_final = 2.0*math.pi/60.0 #rad/s (1 RPM)
    omega_calc = rot_velo[-1]/(60.0*2.0*math.pi)
    I_stage_3_rot = 0.5*m_vehicle*0.15**2.0

    m_yoyo = 1.0 #kg
    number_of_yoyos = 2.0
    arm_yoyo = 1.0*10**-3

    while omega_calc > omega_final:
            
        I_yoyo = number_of_yoyos*m_yoyo*(arm_yoyo**2.0)
        I_total_3 = I_stage_3_rot + I_yoyo

        omega_calc = (((rot_velo[-1]/(60.0*2.0*math.pi))*I_stage_3_rot)/I_total_3)
        arm_yoyo = arm_yoyo + 1.0*10**-3

    print("\nYoyo despin arm length is: {0} mm".format(arm_yoyo*1000.0))

    #Calculate the yoyo mass diameter
    rho_yoyo = 2700.0
    t_yoyo = 20.0/1000.0

    r_yoyo = (((m_yoyo/rho_yoyo)/t_yoyo)/math.pi)**0.5
    d_yoyo = r_yoyo*2.0

    print("Yoyo mass diameter: {0} mm".format(d_yoyo*1000.0))

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