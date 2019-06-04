import math

#Design the spin stabilized controls


def main():

    R_airframe = 150.0/1000.0

    rho_grain_NG_NC = 1133.0
    rho_airframe = 8000.0

    m_propellant = 350.0*0.32
    m_airframe = 100.0*0.32

    I_prop = 0.5*m_airframe*R_airframe**2.0
    I_airframe = m_airframe*R_airframe**2.0

    I_total = I_prop + I_airframe

    v_side = 5.0 #Wind forces in m/s

    rho_altitude = 1.22
    Maximum_downrange_distance = 20000.0 #meters
    h_apogee = 208000.0 #meters

    t_total = 227.0 #seconds
    t_stage = 9.0 #seconds

    C_drag_cylinder = 1.2

    D_core = 300.0/1000.0
    L_core = 270.0/1000.0
    L_fairing = (180.0*6.0)/1000.0

    L_nozzle_1 = 435.0/1000.0
    L_nozzle_2 = 421.0/1000.0
    L_nozzle_3 = 494.0/1000.0

    number_cores_2 = 3.0
    number_cores_3 = 2.0

    L_booster_2 = number_cores_2 * L_core
    L_booster_3 = number_cores_3 * L_core

    A_side = D_core * (L_booster_2 + L_nozzle_2 + L_booster_3 + L_nozzle_3 + L_fairing)

    dt = 1.0/1000.0
    t_flight = 0.0
    h_altitude = 0.0
    a_average = (h_apogee/0.5)/(t_total**2.0)

    L_downrange = 0.0

    omega_p = 0.0 #gyroscopic precession in rad/s
    omega_rot = 2.0*math.pi*100.0 #rad/s

    #Determine the rotation velocity
    while (t_flight < t_total):
        
        h_altitude = 0.5*a_average*t_flight**2.0
        v_altitude = a_average*t_flight

        if (h_altitude < 30000.0):
            rho_altitude, P_altitude, T_altitude = compute_density(h_altitude)

        else:
            rho_altitude = 0.0

        F_side = C_drag_cylinder*A_side*0.5*rho_altitude*v_side**2.0
        omega_p = (F_side * (L_booster_2 + L_nozzle_2 + L_booster_3 + L_nozzle_3 + L_fairing)/2.0)/(I_total * omega_rot)
        total_angle = omega_p*dt

        L_downrange = L_downrange + total_angle*(180.0/math.pi)*v_altitude

        t_flight = t_flight + dt


    print("The total downrange distance is {0} km".format(L_downrange/1000.0))

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