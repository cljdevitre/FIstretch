import unittest
import pandas as pd
import RelaxiFI as relax


# decimalPlace=4
# class test_kmdepth(unittest.TestCase):
#     def test_kmdepth(self):
#         self.assertAlmostEqual(round(relax.find_P_for_kmdepth(target_depth=5, 
# config=relax.config_crustalmodel(crust_dens_kgm3=2750),
# initial_P_guess=0, tolerance=0.1)[0],4), 1.3488,
# decimalPlace, "Calculated P doesnt match test value")
        

decimalPlace=4
class test_SP94_EOS(unittest.TestCase):
    def test_SP94_pressure_to_density(self):
        self.assertAlmostEqual(relax.calculate_rho_for_P_T(P_kbar=1, 
T_K=1400, EOS='SP94')[0], 0.297958,
decimalPlace, "Calculated SP94 P doesnt match test value")
        



