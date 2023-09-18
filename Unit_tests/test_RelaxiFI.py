import unittest
import RelaxiFI as relax

decimalPlace = 4

class TestFindP4kmD(unittest.TestCase):
    def test_find_P_for_kmdepth(self):
        self.assertAlmostEqual(
            round(relax.find_P_for_kmdepth(target_depth=5, 
            config=relax.config_crustalmodel(crust_dens_kgm3=2750),
            initial_P_guess=0, tolerance=0.1)[0],4), 
            float(1.3489),
            decimalPlace,
            "Calculated P doesn't match test value"
        )