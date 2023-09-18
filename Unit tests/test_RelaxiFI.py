import unittest
import pandas as pd
import RelaxiFI as relax


decimalPlace=4
class test_kmdepth(unittest.TestCase):
    def test_kmdepth(self):
        self.assertAlmostEqual(relax.find_P_for_kmdepth(target_depth=5, 
config=relax.config_crustalmodel(crust_dens_kgm3=2750),
initial_P_guess=0, tolerance=0.1)[0], 1.3488,
decimalPlace, "Calculated P doesnt match test value")
        



