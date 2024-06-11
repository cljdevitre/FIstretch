import unittest

def add(a, b):
    return a + b

class TestAddition(unittest.TestCase):

    def test_addition(self):
        result = add(4, 4)
        self.assertEqual(result, 8)

        



