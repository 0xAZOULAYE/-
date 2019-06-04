function [m_payload, frac_payload,frac_structural, frac_propellant, Isp_motor,...
            Cd_z,Cd_r,D_body,D_capsule,Length_total,Length_capsule,Length_body,Omega_r,Omega_theta,...
            m_propellant,m_structural,m_total,m_individual_stage_structural, n_stages,...
            m_capsule, m_stages, m_individual_stage,...
            m_individual_stage_structural, n_stages, t_burn, t_stage, F_thrust, h0, massflow_motor, TWR] = define_constants()

            //This function is meant to define all of the initial values which are fundamental to the operation of the LV

    m_payload = 50.0 //kg
    frac_payload = 10/100
    frac_structural = 20/100
    frac_propellant = 1 - frac_payload - frac_structural

    Isp_motor = 215.0 //seconds

    //Determine the aerodynamic properties of the LV
    Cd_z = 0.03
    Cd_r = 1.0

    D_body = 452.0/1000.0 //m
    D_capsule = 452.0/1000.0 //m

    BC_fraction = 1/4 //Body length to capsule length fraction
    Length_total = D_body*20
    Length_capsule = Length_total * BC_fraction
    Length_body = Length_total - Length_capsule

    //Define the control parameters
    Omega_r = 180 // dps
    Omega_theta = 180 // dps

    //Determine the LV masses
    m_total = 620.0 //1/frac_payload * m_payload
    m_propellant = 370.0 //m_total * frac_propellant
    m_structural = 200.0 //m_total * frac_structural

    //Compute the total burn time as a function of the specific impulse and the required thrust
    h0 = 0 //Sea level altitude

    TWR = 11.0 //Thrust to weight ratio
    F_thrust = [52000.0; 21000.0; 12600.0]
    massflow_motor = [26.0; 5.36; 1.85] //F_thrust / (get_gravity(h0)*Isp_motor)
    t_burn = 12.21 // + 12.09 + 16.77//m_propellant/massflow_motor
    
    //------------------------------------
        
    P_max = 3.7*10.0**6.0
    n_motor = 0.64
    a_motor = (2.3*10**-3)/(1.0*10.0**6)**n_motor
        
    D_max = D_body           
    L_max = 10.0            
    m_max = 397.0           
        
    OD_stage = D_max
    ID_stage = 50.0e-3
        
    r_motor = a_motor*P_max**n_motor
    W_gap = 10.0e-3

    t_stage = [12.21; 12.09; 14.41]//((OD_stage/2 -W_gap)/2.0)/r_motor
        
    n_stages = 1.0//t_burn/t_stage
    
        //-----------------------------------

    m_capsule = m_payload * (1 + frac_structural)
    m_stages = m_total - m_capsule
    m_individual_stage = [m_stages/16*8; m_stages/16*5; m_stages/16*3;]  
    m_individual_stage_structural = [m_individual_stage(1)*frac_structural;m_individual_stage(2)*frac_structural;m_individual_stage(3)*frac_structural]

    disp("The payload mass has been determined at: " + string(m_payload) + " kg")
    disp("The propellant mass has been determined at: " + string(m_propellant) + " kg")
    disp("The structural mass has been determined at: " + string(m_structural) + " kg")
    disp("The total mass has been determined at: " + string(m_total) + " kg" + ascii(10))

    disp("\nThe capsule diameter has been determined at: " + string(D_capsule) + " m")
    disp("The capsule length has been determined at: " + string(Length_capsule) + " m")
    disp("The body diameter has been determined at: " + string(D_body) + " m")
    disp("The body length has been determined at: " + string(Length_body) + " m")
    disp("The total length has been determined at: " + string(Length_total) + " m" + ascii(10))

    disp("\nThe capsule mass has been determined at: " + string(m_capsule) + " kg" + ascii(10))
    disp("The stage mass has been determined at: " + string(m_individual_stage) + " kg" + ascii(10))
    disp("The stage 1 structural mass has been determined at: " + string(m_individual_stage_structural(1)) + " kg" + ascii(10))
    //disp("The stage 2 structural mass has been determined at: " + string(m_individual_stage_structural(2)) + " kg" + ascii(10))
    //disp("The stage 3 structural mass has been determined at: " + string(m_individual_stage_structural(3)) + " kg" + ascii(10))

endfunction


function F_dz = compute_axial_drag()

endfunction


function F_dr = compute_lateral_drag(Length_total, D_body, Cd_r)

    V_wind_max = 63/3.6 //Maximum allowable wind speed in m/s
    Alpha_max = 30/(180*%pi) //Maximum allowable angle during launch

endfunction


function [rho_altitude, P_altitude, T_altitude] = compute_density(altitude)

    T0 = 293.0 //K
    P0 = 101250.0 //Pa
    Rho_0 = 1.225

    R_ideal = 8.31447 // J/mol*K
    M_air = 0.0289644 // kg/mol
    g0 = 9.81 //m/s2

    L = 0.0065 // K/m
    T_altitude = T0 - L*altitude
    P_altitude = P0 * (1.0 - (L*altitude)/T0) ** ((get_gravity(altitude) * M_air) / (R_ideal * L))
    rho_altitude = (P_altitude * M_air) / (R_ideal * T_altitude)

endfunction


//This function is meant to compute the gravitational pull as a function of the altitude
function g1 = get_gravity(altitude)

    G_grav = 6.67408 * 10**-11
    M_earth = 5.972*10**24
    R_earth = 6371.0*10**3

    g1 = (G_grav * M_earth)/((R_earth + altitude)**2)

endfunction

//Obtain the starting values
[m_payload, frac_payload,frac_structural, frac_propellant, Isp_motor,...
Cd_z,Cd_r,D_body,D_capsule,Length_total,Length_capsule,Length_body,Omega_r,Omega_theta,...
                                            m_propellant,m_structural,m_total,m_individual_stage_structural, n_stages,...
                                            m_capsule, m_stages, m_individual_stage,...
                                            m_individual_stage_structural, n_stages, t_burn, t_stage, F_thrust, h0, massflow_motor, TWR] = define_constants()


//Compute within a while loop the velocity and altitude increase within the burn time as a function of a constant thrust
t_flight = 0
m_vehicle = m_total
v_vehicle = 0
a_vehicle = 0
h_vehicle = h0

P_1 = 101250.0
T_1 = 293.0
rho_1 = 1.225

t_stage_computed = 0.0
n_stages_burnt = 1.0

//Define the time step
dt = 1/1000

F_1 = 52402.0
Isp_motor = 212.0
m_1 = F_1/(Isp_motor*9.81)

F_thrust = [F_1; 21000.0; 16000.0]
massflow_motor = [30.0; 7.24; 4.6] //F_thrust / (get_gravity(h0)*Isp_motor)
t_stage = [12.21; 12.09; 14.41]//((OD_stage/2 -W_gap)/2.0)/r_motor

counter_thrust = 1


//Compute the staging of the rocket
while n_stages_burnt <= (n_stages)

    t_stage_computed = 0.0

    disp("Stage " + string(n_stages_burnt) + " thrust is equal to " + string(F_thrust(counter_thrust)) + " N")

    while t_stage_computed < t_stage(n_stages_burnt)

        //Compute the current vehicle mass
        T0 = 293.0 //K
        P0 = 101250.0 //Pa
        Rho_0 = 1.225
    
        R_ideal = 8.31447 // J/mol*K
        M_air = 0.0289644 // kg/mol
        g0 = 9.81 //m/s2
    
        L = 0.0065 // K/m
        T_altitude = T0 - L*h_vehicle
        P_altitude = P0 * (1.0 - (L*h_vehicle)/T0) ** ((get_gravity(h_vehicle) * M_air) / (R_ideal * L))
        rho_altitude = (P_altitude * M_air) / (R_ideal * T_altitude)
        
        a_vehicle = (F_thrust(counter_thrust) - get_gravity(h_vehicle)*m_vehicle - 0.5*(D_capsule/2)**2*%pi*Cd_z*rho_altitude*v_vehicle**2)/m_vehicle

        m_vehicle = m_vehicle - massflow_motor(counter_thrust)*dt
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        t_flight = t_flight + dt
        t_stage_computed = t_stage_computed + dt

    end
    

    m_vehicle = m_vehicle - m_individual_stage_structural(counter_thrust)

    n_stages_burnt = n_stages_burnt + 1.0
    counter_thrust = counter_thrust + 1

    disp(string(n_stages_burnt))

end

//disp("The total burn time is equal to " + string(t_burn) + " seconds")
disp("The burnout density is equal to " + string(rho_altitude) + " kg/m3" + ascii(10))
disp("The burnout pressure is equal to " + string(P_altitude) + " Pa" + ascii(10))
disp("The burnout temperature is equal to " + string(T_altitude) + " K\n" + ascii(10))

disp("The stage burn time is equal to " + string(t_stage) + " seconds")

disp("The burnout altitude is equal to " + string(h_vehicle) + " meters")
disp("The burnout velocity is equal to " + string(v_vehicle) + " m/s" + ascii(10))


while v_vehicle > 0.0

        //Compute the current vehicle mass
                //Compute the current vehicle mass
        T0 = 293.0 //K
        P0 = 101250.0 //Pa
        Rho_0 = 1.225
    
        R_ideal = 8.31447 // J/mol*K
        M_air = 0.0289644 // kg/mol
        g0 = 9.81 //m/s2
    
        if rho_altitude > 0.000000001

            L = 0.0065 // K/m
            T_altitude = T0 - L*h_vehicle
            P_altitude = P0 * (1.0 - (L*h_vehicle)/T0) ** ((get_gravity(h_vehicle) * M_air) / (R_ideal * L))
            rho_altitude = (P_altitude * M_air) / (R_ideal * T_altitude)

        end

        a_vehicle = (0.0 - get_gravity(h_vehicle)*m_vehicle - 0.5*(D_capsule/2)**2*%pi*Cd_z*rho_altitude*v_vehicle**2)/m_vehicle
    
        m_vehicle = m_vehicle - 0.0
        v_vehicle = v_vehicle + a_vehicle*dt
        h_vehicle = h_vehicle + v_vehicle*dt

        t_flight = t_flight + dt

end

disp("The apogee altitude is equal to " + string(h_vehicle) + " meters")
disp("The apogee velocity is equal to " + string(v_vehicle) + " m/s")
disp("The apogee time is equal to " + string(t_flight) + " s")

//Determine the distance to the horizon
r_akkal = 6371e3 //m
r_timezdith = r_akkal + h_vehicle

phi_akkal= acos(r_akkal/r_timezdith)
l_igli = phi_akkal*r_akkal

disp("The horizon distance is:" + string(l_igli/1000.0) + " km")

exit