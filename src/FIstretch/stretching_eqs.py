import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from scipy.optimize import newton
import warnings

from FIstretch.density_depth_crustal_profiles import *
from FIstretch.EOS_CO2 import *

## Functions to find P when the user chooses to start with a depth. It requires input of a crustal model

class config_crustalmodel:
    """
    A configuration class for specifying parameters of the crustal model. 

    Attributes:
    - crust_dens_kgm3 (float): The density of the crust in kilograms per cubic meter (kg/m^3).
    - d1 (float): The depth boundary for the first layer in kilometers (km).
    - d2 (float): The depth boundary for the second layer in kilometers (km).
    - rho1 (float): The density of the first layer in kilograms per cubic meter (kg/m^3).
    - rho2 (float): The density of the second layer in kilograms per cubic meter (kg/m^3).
    - rho3 (float): The density of the third layer in kilograms per cubic meter (kg/m^3).
    - model (str): The name of the model used for crustal calculations.
    """
    def __init__(self, crust_dens_kgm3=None,
                 d1=None, d2=None, rho1=None, rho2=None, rho3=None, model=None):
        self.crust_dens_kgm3 = crust_dens_kgm3
        self.d1 = d1
        self.d2 = d2
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3
        self.model = model

def objective_function_depth(P_kbar, target_depth, crust_dens_kgm3,
                            d1, d2, rho1, rho2, rho3, model):
    """
    Calculate the difference between the current depth and the target depth
    given pressure (P_kbar) and other parameters.

    Parameters:
    - P_kbar (float): The pressure in kilobars (kbar) to be used in the depth calculation.
    - target_depth (float): The desired depth in kilometers (km).
    - crust_dens_kgm3 (float): The density of the crust in kilograms per cubic meter (kg/m^3).
    - d1, d2 (float): Depth boundaries for different layers (km).
    - rho1, rho2, rho3 (float): Densities for different layers (kg/m^3).
    - model (str): The name of the model used for the depth calculation.

    Returns:
    - float: The difference between the current depth and the target depth.
    """

    current_depth = convert_pressure_to_depth(P_kbar=P_kbar, crust_dens_kgm3=crust_dens_kgm3, g=9.81,
                                              d1=d1, d2=d2, rho1=rho1, rho2=rho2, rho3=rho3, model=model)[0]

    return current_depth - target_depth

def find_P_for_kmdepth(target_depth, config=config_crustalmodel(), initial_P_guess=0, tolerance=0.1):
    """
    Approximate the pressure (P_kbar) based on the target depth using the Newton-Raphson method.

    Parameters:
    - target_depth (float, Pandas Series, list): The desired depth(s) in kilometers (km).
    - P_kbar (float, optional): Initial guess for the pressure in kilobars (kbar). Default is None.
    - crust_dens_kgm3 (float, optional): The density of the crust in kilograms per cubic meter (kg/m^3). Default is None.
    - d1, d2 (float, optional): Depth boundaries for different layers (km). Default is None.
    - rho1, rho2, rho3 (float, optional): Densities for different layers (kg/m^3). Default is None.
    - model (str, optional): The name of the model used for the depth calculation. Default is None.
    - tolerance (float, optional): How close the pressure estimate should be to the true value. Default is 0.1.

    Returns:
    - float or Pandas Series or list: The estimated pressure(s) (P_kbar) that correspond to the target depth(s).
    """

    if isinstance(target_depth, (float, int)):
        target_depth = [target_depth]  

    pressures = []

    for depth in target_depth:
        if all(v is None for v in [config.crust_dens_kgm3, config.d1, config.d2, config.rho1, config.rho2, config.rho3, config.model]):
            config.crust_dens_kgm3 = 2750
            warning_message = "\033[91mNo crustal parameters were provided, setting crust_dens_kgm3 to 2750. \nPlease use config_crustalmodel(...) to set your desired crustal model parameters.\033[0m"
            warnings.simplefilter("always")
            warnings.warn(warning_message, Warning, stacklevel=2)
        
        # Use the Newton-Raphson method for each target depth
        pressure = newton(objective_function_depth, initial_P_guess, args=(depth, config.crust_dens_kgm3, config.d1, config.d2, config.rho1, config.rho2, config.rho3, config.model), tol=tolerance)
        pressures.append(pressure)

    if isinstance(target_depth, (float, int)):
        return pressures[0]
    elif isinstance(target_depth, pd.Series):
        return pd.Series(pressures)
    else:
        return pressures

## Functions to calculate P and rho from ideal gas law

# Density and mass
def ideal_calc_rho_for_r_P_T(*,P_MPa,T,r):
    """
    Calculate the density and mass of CO2 using the ideal gas law for a sphere of radius r.

    Parameters:
    P_MPa (float): Pressure in MegaPascals (MPa).
    T (float): Temperature in Kelvin (K).
    r (float): Radius of the sphere in meters (m).

    Returns:
    tuple: A tuple containing two values.
        - rho (float): Density of the gas in grams per cubic centimeter (g/cm^3).
        - m (float): Mass of CO2 in kilograms (kg).
    """
    P=P_MPa*10**6 # convert MPa to Pa
    M=44.01/1000 #kg/mol
    V=4/3*math.pi*r**3 #m3
    R=8.314 #J.mol/K J: kg·m²/s² 
    m=P*V*M/(R*T) #CO2 mass in kg
    rho=(m/V)/1000 # rho in g/cm3
    return rho, m/1000

# Pressure
def ideal_calc_P_for_V_rho_T(*,co2_mass_g,T,r):
    """
    Calculate the pressure of CO2 using the ideal gas law for a sphere of radius r.

    Parameters:
    co2_mass_g (float): Mass of CO2 gas in grams (g).
    T (float): Temperature in Kelvin (K).
    r (float): Radius of the sphere in meters (m).

    Returns:
    tuple: A tuple containing two values.
        - rho (float): Density of the gas in grams per cubic centimeter (g/cm^3).
        - P_MPa (float): Pressure in MegaPascals (MPa).

    """
    M=44.01/1000 #kg/mol
    m=co2_mass_g*1000
    V=4/3*math.pi*r**3 #m3
    R=8.314 #J.mol/K J: kg·m²/s² 
    P=m*R*T/(M*V) #P in Pa
    P_MPa=P/(10**6) #P in MPa
    rho=(m/V)/1000 #rho in g/cm3
    return rho,P_MPa

## Functions to calculate initial conditions

# Initial volume, CO2 density and CO2 mass in FI
def calculate_initial_V_CO2rho_mass(*,EOS='SW96',P,T,r):
    """
    Calculate the initial volume, CO2 density, and CO2 mass inside a fluid inclusion.

    Parameters:
    - EOS (str): Equation of State to use for calculations. Default is 'SW96'.
    - P (float): Pressure inside the fluid inclusion in MPa.
    - T (float): Temperature inside the fluid inclusion in Kelvin.
    - r (float): Radius of the fluid inclusion in meters.

    Returns:
    - V (float): Volume of the fluid inclusion in cm^3 (assumes a spherical shape).
    - CO2_dens (float): CO2 density inside the fluid inclusion in g/cm^3.
    - CO2_mass (float): CO2 mass inside the fluid inclusion in grams.
    """
    r=r*10**2
    V=4/3*math.pi*r**3 #cm3, Volume of the FI, assume sphere
    P_kbar=P/100 #Internal pressure of the FI
    CO2_dens=calculate_rho_for_P_T(EOS=EOS,P_kbar=P_kbar,T_K=T)[0] #g/cm3, CO2 density, calc by Span&Wagner(96)
    CO2_mass=CO2_dens*V # this is our CO2 mass in the FI
    return V, CO2_dens, CO2_mass


def calculate_DPdt(ascent_rate_ms,config=config_crustalmodel(),D_initial=None,D_final=None,D_step=100,initial_P_guess=0, tolerance=0.001):
    """
    Calculate the decompression rate (DP/dt) during ascent.

    Parameters:
    - ascent_rate_ms (float): Ascent rate in meters per second.
    - D_initial (float): Initial depth in kilometers. Default is 30 km.
    - D_final (float): Final depth in kilometers. Default is 0 km.
    - D_step (int): Number of depth steps for calculation. Default is 100.

    Returns:
    - D (pd.Series): Depth values in kilometers.
    - Pexternal_steps (list): Lithostatic pressure values in MPa at each depth step.
    - dt (float): Time step for the integration.
    """

    if D_initial is None or D_final is None or D_initial <= D_final:
        raise ValueError("Both D_initial and D_final must be provided, and D_initial must be larger than D_final")
    if D_initial>30 and D_step <= 80 and ascent_rate_ms <= 0.02:
        raise Warning("Your D_step is too small, the minimum recommended for ascent rates below 0.02 m/s is 80")
    D = pd.Series(list(np.linspace(D_initial, D_final, D_step)))  # km

    Pexternal_steps=find_P_for_kmdepth(D, config=config, initial_P_guess=initial_P_guess, tolerance=tolerance)
    Pexternal_steps_MPa=Pexternal_steps*100

    # Time steps of the ascent
    ascent_rate = ascent_rate_ms / 1000  # km/s
    D_change = abs(D.diff())
    time_series = D_change / ascent_rate  # calculates the time in between each step based on ascent rate
    dt = time_series.max()  # this sets the time step for the iterations later

    return D, Pexternal_steps_MPa, dt


## Auxilliary functions for the stretching model

# Helper function for internal pressure and co2 density
def calculate_step_P_for_m_r(*,EOS='SW96',m,T,r):
    """
    Calculate the internal pressure, volume, and CO2 density of an inclusion as a function of temperature (T),
    inclusion radius (r in cm), and CO2 mass (m in g).

    Parameters:
    - EOS (str): Equation of State for CO2 density calculation (default: 'SW96').
    - m (float): CO2 mass in grams.
    - T (float): Temperature in Kelvin.
    - r (float): Radius of the inclusion in meters

    Returns:
    - V (float): Volume of the inclusion in cm^3 (assuming a sphere).
    - CO2_dens (float): CO2 density in g/cm^3.
    - P_new (float): Internal pressure in MPa.
    """
    r=r*10**2
    V=4/3*math.pi*r**3 #cm3, Volume of the FI, assume sphere
    CO2_dens=m/V
    try:
        P_new=calculate_P_for_rho_T(EOS=EOS,CO2_dens_gcm3=CO2_dens, T_K=T)['P_MPa'][0] #g/cm3, CO2 density, calc by Span&Wagner(96)
        return V, CO2_dens, P_new
    except ValueError:
        return V,CO2_dens,np.nan

class power_creep_law_constants:
    """
    Olivine power-law creep constants used in the stretching model (Wanamaker and Evans, 1989).

    Attributes:
    - A (float): Creep law constant A (default: 3.9e3).
    - n (float): Creep law constant n (default: 3.6).
    - Q (float): Activation energy for dislocation motions in J/mol (default: 523000).
    - IgasR (float): Gas constant in J/(mol*K) (default: 8.314).
    """
    def __init__(self):
        self.A = 3.9*10**3 #7.0 * 10**4
        self.n = 3.6 #3
        self.Q = 523000 # 520 Activation energy for dislocation motions in J/mol
        self.IgasR= 8.314  # Gas constant in J/(mol*K)

# Helper function to calculate dR/dt
def calculate_dR_dt(*,R, b, T,  Pinternal, Pexternal):
    """
    Calculate the rate of change of inclusion radius (dR/dt) based on power law creep.

    Parameters:
    - R (float): Inclusion radius in m.
    - b (float): Distance to the crystal defect structures, Wanamaker and Evans (1989) use R0/b=1/1000
    - T (float): Temperature in Kelvin.
    - Pinternal (float): Internal pressure in MPa.
    - Pexternal (float): External pressure in MPa.

    Returns:
    - dR_dt (float): Rate of change of inclusion radius in m/s.
    """
    pl_Cs = power_creep_law_constants()
    if Pinternal<Pexternal==True:
        S=-1
    else:
        S=1
    try:
        dR_dt = 2 * (S * pl_Cs.A * math.exp(-pl_Cs.Q / (pl_Cs.IgasR * T))) * (((R * b)**3) / (((b**(3 / pl_Cs.n)) - (R**(3 / pl_Cs.n))))**pl_Cs.n) * (((3 * abs(Pinternal - Pexternal)) / (2 * pl_Cs.n))**pl_Cs.n) / R**2
        return dR_dt

    except FloatingPointError:
        return np.nan
    

## Stretching Models

# This function is to model FI stretching during decompression and ascent 
def stretch_in_ascent(*,R0, b, T,ascent_rate_ms,crustal_model_config=config_crustalmodel(),depth_path_ini_fin_step=[100,0,100],EOS,plotfig=True,display_df=True,initial_P_guess=0,tolerance=0.001):

    D, Pexternal_steps, dt = calculate_DPdt(ascent_rate_ms=ascent_rate_ms,config=crustal_model_config,
                                            D_initial=depth_path_ini_fin_step[0],D_final=depth_path_ini_fin_step[1],
                                            D_step=depth_path_ini_fin_step[2],
                                            initial_P_guess=initial_P_guess, tolerance=tolerance)
    Pinternal = Pexternal_steps[0]
    _,CO2_dens_initial,CO2_mass_initial=calculate_initial_V_CO2rho_mass(EOS=EOS,P=Pinternal,T=T,r=R0)
    R_values = [R0]  # List to store R values at different time points
    Pinternal_list=[Pinternal]
    CO2_dens_list=[CO2_dens_initial]
    dR_dt_list=[]
    t0=0
    t_list=[t0]
    results=pd.DataFrame(columns={'Pexternal(MPa)','Pinternal(MPa)',
                                  'Depth(km)','Fi_radius(\u03BCm)',
                                  'CO2_dens_gcm3','dR/dt(m/s)' })
    for i in range(len(Pexternal_steps)):
        Pexternal = Pexternal_steps[i]
        dR_dt = calculate_dR_dt(R=R_values[-1], b=b,Pinternal=Pinternal, Pexternal=Pexternal, T=T)
        R_new = R_values[-1] + dR_dt * dt
        _,CO2_dens_new,P_new=calculate_step_P_for_m_r(EOS=EOS,m=CO2_mass_initial,T=T,r=R_new) # Why am I multiflying here??
        Pinternal=P_new
        dR_dt_list.append(dR_dt)
        R_values.append(R_new)
        Pinternal_list.append(Pinternal)
        CO2_dens_list.append(CO2_dens_new)
        ti=t_list[-1]+dt
        t_list.append(ti)
    
    results['Pexternal(MPa)']=Pexternal_steps
    results['Pinternal(MPa)']=Pinternal_list[1:]
    results['Depth(km)']=D
    results['Fi_radius(\u03BCm)']=[num * 10**6 for num in R_values][1:]
    results['CO2_dens_gcm3']=CO2_dens_list[1:]
    results['dR/dt(m/s)']=dR_dt_list
    results['\u0394R/R0']=(results['Fi_radius(μm)']-results['Fi_radius(μm)'][0])/results['Fi_radius(μm)'][0]
    results['Time (s)']=[t - dt for t in t_list[1:]]
    if display_df==True:
        display(results.head())
    if plotfig==True:
        fig, (ax0,ax1) = plt.subplots(1,2, figsize=(10,3))
        ax0.plot(-results['Depth(km)'],results['ΔR/R0'],marker='s')
        ax0.set_xlabel("Depth")
        ax0.set_ylabel("DeltaR/R0")

        ax1.plot(-results['Depth(km)'],results['CO2_dens_gcm3'],marker='s')
        ax1.set_xlabel("Depth")
        ax1.set_ylabel("CO2_density_gmL")

    return results, fig if 'fig' in locals() else results

# This function is to model stretching at fixed External Pressure (e.g., during stalling or upon eruption)
def stretch_at_constant_Pext(*,R,b=None,T,EOS='SW96',Pinternal,Pexternal,totaltime,steps,calc_method='Euler',RKorder=4,display_df=True,plotfig=False):

    dt=totaltime/steps

    if b is None:
        b=R*1000
    
    if EOS=='ideal':
        CO2_dens_initial,CO2_mass_initial=ideal_calc_rho_for_r_P_T(P_MPa=Pinternal,T=T,r=R)
    else:
        _,CO2_dens_initial,CO2_mass_initial=calculate_initial_V_CO2rho_mass(EOS=EOS,P=Pinternal,T=T,r=R)

    results = pd.DataFrame([{'Time(s)': 0,
                             'Step':0,
                             'dt(s)':0,
                            'Pexternal(MPa)': Pexternal,
                            'Pinternal(MPa)': Pinternal,
                            'dR/dt(m/s)': calculate_dR_dt(R=R, b=b, Pinternal=Pinternal, Pexternal=Pexternal, T=T),
                            'Fi_radius(μm)': R*10**6,
                            'b (distance to xtal rim -μm)':b*10**6,
                            '\u0394R/R0 (fractional change in radius)':0,
                            'CO2_dens_gcm3': CO2_dens_initial}], index=range(steps))

    if calc_method=='Euler':
        for step in range(1,steps):
            dR_dt = calculate_dR_dt(R=R, b=b,Pinternal=Pinternal, Pexternal=Pexternal, T=T)
            R_new= R + dR_dt*dt
            if EOS=='ideal':
                CO2_dens_new,P_new=ideal_calc_P_for_V_rho_T(co2_mass_g=CO2_mass_initial,T=T,r=R_new)
            else:
                _,CO2_dens_new,P_new=calculate_step_P_for_m_r(EOS=EOS,m=CO2_mass_initial,T=T,r=R_new)
            R=R_new
            b=1000*R
            Pinternal=P_new
            
            results.loc[step] = [step * dt, step, dt, Pexternal, Pinternal, dR_dt, R * 10 ** 6,
                                (R * 10 ** 6 - results.loc[0, 'Fi_radius(μm)']) / results.loc[0, 'Fi_radius(μm)'],
                                b * 10 ** 6, CO2_dens_new]

    if calc_method=='RK':

        for step in range(steps):
            if RKorder == 2:
                k1 = dt * calculate_dR_dt(R=R, b=b, T=T, Pinternal=Pinternal, Pexternal=Pexternal)
                k2 = dt * calculate_dR_dt(R=R + 0.5 * k1, b=b, T=T, Pinternal=Pinternal, Pexternal=Pexternal)

                dR_dt = ((k1 + k2) / 2) / dt
                R += (k1 + k2) / 2
            elif RKorder == 4:
                k1 = dt * calculate_dR_dt(R=R, b=b, T=T, Pinternal=Pinternal, Pexternal=Pexternal)
                k2 = dt * calculate_dR_dt(R=R + 0.5 * k1, b=b, T=T, Pinternal=Pinternal, Pexternal=Pexternal)
                k3 = dt * calculate_dR_dt(R=R + 0.5 * k2, b=b, T=T, Pinternal=Pinternal, Pexternal=Pexternal)
                k4 = dt * calculate_dR_dt(R=R + k3, b=b, T=T, Pinternal=Pinternal, Pexternal=Pexternal)

                dR_dt = ((k1 + 2 * k2 + 2 * k3 + k4) / 6) / dt
                R += (k1 + 2 * k2 + 2 * k3 + k4) / 6

            if EOS == 'ideal':
                CO2_dens_temp, P_temp = ideal_calc_P_for_V_rho_T(co2_mass_g=CO2_mass_initial, T=T, r=R)
            else:
                _, CO2_dens_temp, P_temp = calculate_step_P_for_m_r(EOS=EOS, m=CO2_mass_initial, T=T, r=R)

            Pinternal = P_temp
            b = R * 1000
            results.loc[step] = [step * dt, step, dt, Pexternal, Pinternal, dR_dt, R * 10 ** 6,
                                (R * 10 ** 6 - results.loc[0, 'Fi_radius(μm)']) / results.loc[0, 'Fi_radius(μm)'],
                                b * 10 ** 6, CO2_dens_temp]
    if display_df==True:
        display(results.head())

    if plotfig==True:
        fig, (ax0,ax1) = plt.subplots(1,2, figsize=(10,3))
        ax0.plot(-results['Time(s)'],results['\u0394R/R0 (fractional change in radius)'],marker='s')
        ax0.set_xlabel("Time (s)")
        ax0.set_ylabel("\u0394R/R0 (fractional change in radius)")

        ax1.plot(-results['Time(s)'],results['CO2_dens_gcm3'],marker='s')
        ax1.set_xlabel("Time(s)")
        ax1.set_ylabel("CO2_density_gmL")

    return results, fig if 'fig' in locals() else results


def calculate_results(R_values, b_values, T, EOS, Pinternal, Pexternal, totaltime, steps, T4endcalc_PD, calc_method='Euler',
                      display_df=False, plotfig=False, config=config_crustalmodel(crust_dens_kgm3=2750)):
    results_dict = {}
    print(config.crust_dens_kgm3)
    for idx_R, R in enumerate(R_values):  # Use enumerate to get the index of R_values
        R_key = f'R{idx_R}'  # Use 'R' followed by the index
        results_dict[R_key] = {}

        for idx_b, b in enumerate(b_values):  # Use enumerate to get the index of b_values
            b_key = f'b{idx_b}'  # Use 'b' followed by the index
            results = stretch_at_constant_Pext(R=R, b=b, T=T, Pinternal=Pinternal, Pexternal=Pexternal,
                                              totaltime=totaltime, steps=steps, EOS=EOS, calc_method=calc_method,
                                              display_df=display_df, plotfig=plotfig)
            results = results[0]
            results['Calculated depths (km)_StorageT'] = convert_pressure_to_depth(
                P_kbar=results['Pinternal(MPa)'] / 100,
                crust_dens_kgm3=config.crust_dens_kgm3, g=9.81,
                d1=config.d1, d2=config.d2,
                rho1=config.rho1, rho2=config.rho2, rho3=config.rho3,
                model=config.model)
            results['Calculated P from rho (MPa)_TrappingT'] = calculate_P_for_rho_T(
                EOS='SW96', CO2_dens_gcm3=results['CO2_dens_gcm3'], T_K=T4endcalc_PD + 273.15)['P_MPa']

            results['Calculated depths (km)_TrappingT'] = convert_pressure_to_depth(
                P_kbar=results['Calculated P from rho (MPa)_TrappingT'] / 100,
                crust_dens_kgm3=config.crust_dens_kgm3, g=9.81,
                d1=config.d1, d2=config.d2,
                rho1=config.rho1, rho2=config.rho2, rho3=config.rho3,
                model=config.model)

            results_dict[R_key][b_key] = results

    return results_dict
