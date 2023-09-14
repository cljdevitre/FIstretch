import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from scipy.optimize import newton
import warnings

from FIstretch.density_depth_crustal_profiles import *
from FIstretch.EOS_CO2 import *

###### This set of functions is used to find P when the user chooses to start with a depth. It requires input of a crustal model
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

def find_P_for_kmdepth(target_depth, config=config_crustalmodel(), initial_P_guess=0,tolerance=0.1):
    """
    Approximate the pressure (P_kbar) based on the target depth using the Newton-Raphson method.

    Parameters:
    - target_depth (float): The desired depth in kilometers (km).
    - P_kbar (float, optional): Initial guess for the pressure in kilobars (kbar). Default is None.
    - crust_dens_kgm3 (float, optional): The density of the crust in kilograms per cubic meter (kg/m^3). Default is None.
    - d1, d2 (float, optional): Depth boundaries for different layers (km). Default is None.
    - rho1, rho2, rho3 (float, optional): Densities for different layers (kg/m^3). Default is None.
    - model (str, optional): The name of the model used for the depth calculation. Default is None.
    - tolerance (float, optional): How close the pressure estimate should be to the true value. Default is 0.1.

    Returns:
    - float: The estimated pressure (P_kbar) that corresponds to the target depth.
    """

    if all(v is None for v in [config.crust_dens_kgm3, config.d1, config.d2, config.rho1, config.rho2, config.rho3, config.model]):
        config.crust_dens_kgm3 = 2750
        warnings.warn("No crustal parameters were provided, setting crust_dens_kgm3 to 2750. Please use config_crustalmodel(...) to set your desired crustal model parameters.", UserWarning)

    pressure = newton(objective_function_depth, initial_P_guess, args=(target_depth, config.crust_dens_kgm3, config.d1, config.d2, config.rho1, config.rho2, config.rho3, config.model), tol=tolerance)
    return pressure