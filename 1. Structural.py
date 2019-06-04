import math


def main():

    D_airframe = 452.0/1000.0
    L_airframe = 6.0

    D_max = 432.0/1000.0
    Cd_side = 1.143

    mat_airframe_name = "304"
    mat_airframe_E = 210.0e9
    mat_airframe_sigma_ys = 215.0e6
    mat_airframe_sigma_shear = mat_airframe_sigma_ys*0.7
    mat_airframe_sigma_fat = 140.0e6
    mat_airframe_max_T = 755.0 #K
    mat_airframe_CT = 1000.0 #Thermal capacity
    mat_airframe_rho = 8000.0
    mat_airframe_sigma_design = 83.0e6 #Design stress at 482 celsius

    mat_mech_name = "7075-T6"
    mat_mech_E = 69.0e9
    mat_mech_sigma_ys = 503.0e6
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
    mat_insulation_rho = 2580.0
    mat_insulation_k = 0.05 #W/m.k

    mat_fairing_name = "ABS"
    mat_fairing_E = 2.6e9
    mat_fairing_sigma_ys = 113.0e6
    mat_fairing_sigma_shear = mat_fairing_sigma_ys*0.6
    mat_fairing_sigma_fat = 0.0
    mat_fairing_max_T = 393.0 #K
    mat_fairing_CT = 0.0 #Thermal capacity
    mat_fairing_rho = 1200.0 #Thermal capacity


    #Define the airframe loadings
    F_thrust_max = 52.0*10**3
    P_max = 5500000.0
    L_cm = 2.89
    L_cp = 2.35

    m_max = 500.0
    g0 = 9.81

    L_fin = 0.30
    t_fin = 5.0e-3
    omega_airframe = 160.0 #RPM
    alpha_airframe = 1.0/(180.0/math.pi) #degrees

    sigma_vm = mat_mech_sigma_shear*2.0
    t_wall = 1.0/1000.0

    sigma_11 = 0.0
    sigma_22 = 0.0
    sigma_33 = 0.0
    
    sigma_12 = 0.0
    sigma_23 = 0.0
    sigma_13 = 0.0
    
    sigma_bend = 0.0

    while sigma_vm > mat_mech_sigma_ys:

        #Compute the tensile stresses
        F_pressure = math.pi*(D_airframe/2.0)**2*P_max
        F_11 = F_thrust_max

        #Compute the bending stresses
        I_x = (math.pi/4.0)*((D_airframe/2.0 + t_wall)**4.0 - (D_airframe/2.0)**4.0)

        M_grav = m_max*g0*(L_cm - L_cp)* math.sin(alpha_airframe)
        M_q = (L_cm - L_cp) * 0.35 * 8.0 * Cd_side * 800000.0 * math.sin(alpha_airframe)
        
        sigma_bend = ((M_grav + M_q)*D_airframe/2.0)/I_x

        sigma_11 = (F_11 + F_pressure)/(math.pi*(D_airframe/2.0 + t_wall)**2.0 - math.pi*(D_airframe/2.0)**2.0) + sigma_bend
        sigma_22 = (P_max*(D_airframe/2.0))/t_wall

        #Compute the shear torque from the despin
        omega = (omega_airframe/60.0*(360.0))/(180.0/math.pi)
        t_despin = 0.1
        alpha_spin = omega/t_despin

        I_z = 0.5*m_max*(D_max/2.0)**2.0
        T_despin = I_z*alpha_spin

        J_torsion = math.pi * ((D_airframe + t_wall)**4.0 - D_airframe**4.0)/32.0
        Sigma_12 = (T_despin * D_airframe/2.0)/J_torsion

        sigma_vm = (0.5*((sigma_11 - sigma_22)**2 + (sigma_22 - sigma_33)**2 + (sigma_33 - sigma_11)**2 + 6.0*(sigma_12**2 + sigma_23**2 + sigma_13**2)))**0.5

        t_wall = t_wall + 0.5/10000.0


    t_horizontal_plate = t_wall


    #Compute the nozzle thickness    
    sigma_vm_nozzle = mat_airframe_sigma_fat*2.0
    t_wall_nozzle = 1.0/1000.0
    D_nozzle = 420.0/1000.0
    L_nozzle = 800.0/1000.0

    sigma_11_nozzle = 0.0
    sigma_22_nozzle = 0.0
    sigma_33_nozzle = 0.0
    
    sigma_12_nozzle = 0.0
    sigma_23_nozzle = 0.0
    sigma_13_nozzle = 0.0

    while sigma_vm_nozzle > mat_airframe_sigma_ys:

        F_pressure_nozzle = math.pi*(D_nozzle/2.0)**2*P_max
        F_11_nozzle = F_thrust_max

        sigma_11_nozzle = (F_11_nozzle + F_pressure_nozzle)/(math.pi*(D_nozzle/2.0 + t_wall_nozzle)**2.0 - math.pi*(D_nozzle/2.0)**2.0)
        sigma_22_nozzle = (P_max*(D_nozzle/2.0))/t_wall_nozzle

        sigma_vm_nozzle = (0.5*((sigma_11_nozzle - sigma_22_nozzle)**2 + (sigma_22_nozzle - sigma_33_nozzle)**2 + (sigma_33_nozzle - sigma_11_nozzle)**2 + 6.0*(sigma_12_nozzle**2 + sigma_23_nozzle**2 + sigma_13_nozzle**2)))**0.5

        t_wall_nozzle = t_wall_nozzle + 0.5/1000.0

    #Check the airframe for buckling

    #Compute the mass of the airframe

    V_airframe = (math.pi*(D_airframe/2.0 + t_wall)**2.0-math.pi*(D_airframe/2.0)**2.0)*L_airframe + (math.pi*(D_max/2.0)**2.0 - math.pi*(D_airframe/2.0)**2.0)*t_horizontal_plate*2.0*(2.0*3.0)
    m_nozzle = (math.pi*(D_nozzle/2.0 + t_wall_nozzle)**2.0 - math.pi*(D_nozzle/2.0)**2.0)*L_nozzle*8000.0
    m_airframe = V_airframe*mat_mech_rho + m_nozzle*3.0

    D_rivet = 9.5/1000.0
    sigma_mat_rivet_ys = 165000000.0
    sigma_mat_rivet_shear = sigma_mat_rivet_ys * 0.8

    number_of_rivets_per_side = 1.0
    sigma_shear_rivets = sigma_mat_rivet_shear * 2.0

    while sigma_shear_rivets > sigma_mat_rivet_shear:
            
        F_hoop = P_max * 0.29 * D_airframe * 0.5
        sigma_shear_rivets = (F_hoop/number_of_rivets_per_side) / (math.pi * (D_rivet/2.0)**2.0)

        number_of_rivets_per_side = number_of_rivets_per_side + 1.0

    print("Computed number of rivets per side: {0} rivets".format(number_of_rivets_per_side))

    number_of_rivets_per_side = 1.0
    sigma_shear_rivets = sigma_mat_rivet_shear * 2.0

    while sigma_shear_rivets > sigma_mat_rivet_shear:
            
        F_z = P_max * math.pi* (D_airframe/2.0)**2.0
        sigma_shear_rivets = (F_z/number_of_rivets_per_side) / (math.pi * (D_rivet/2.0)**2.0)

        number_of_rivets_per_side = number_of_rivets_per_side + 1.0

    #Compute the interstage structure
    number_of_rods = 8.0
    alpha_rod = 9.0/(180.0/math.pi)
    k_factor = 1.2
    l_rod = 0.6
    R_rod_internal = 5.0*10**-3

    F_linear = F_thrust_max*math.cos(alpha_rod)/(number_of_rods)

    I_rod_min = (F_linear*(k_factor*l_rod)**2.0)/((math.pi**2.0)*mat_airframe_E)
    R_outer_rod = (I_rod_min/(math.pi/4.0) + R_rod_internal**4.0)**(1.0/4.0)

    print("Computed interstage buckling rod outer diameter: {0} mm".format(R_outer_rod*2000.0))

    A_rod = F_linear/(mat_mech_sigma_ys)
    R_outer_rod = ((A_rod + math.pi*R_rod_internal**2.0)/math.pi)**0.5

    print("Computed interstage tensile rod outer diameter: {0} mm".format(R_outer_rod*2000.0))


    print("Computed number of rivets per top: {0} rivets".format(number_of_rivets_per_side))

    print("Computed airframe wall thickness: {0} mm".format(t_wall*1000.0))
    print("Computed nozzle wall thickness: {0} mm".format(t_wall_nozzle*1000.0))


    print("Computed airframe mass: {0} kg".format(m_airframe))

main()
