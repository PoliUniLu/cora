import unittest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import cora
import pandas

class KnownValues(unittest.TestCase):
    known_values=[
        [
            [
                ('K', {0,1}),
                ('L', {0,2}),
                ('M', {1,5}),
                ('N', {2,6}),
                ('P', {5,7}),
                ('Q', {6,7})
            ],
            {0,1,2,5,6,7},
            [['K', 'L', 'P', 'Q'],
             ['K', 'M', 'N', 'Q'],
             ['K', 'N', 'P'],
             ['L', 'M', 'N', 'P'],
             ['L', 'M', 'Q']] 
        ]

    ]
    
   

    
    def test_petric_function(self):
        for implicants, coverage,irr_sums in self.known_values:
            result=cora._find_irredundant_sums(implicants, coverage)
            self.assertEqual(irr_sums,result)


class KnownValues(unittest.TestCase):
    known_values = [
        [
        [("P1", {8}),
          ("P2", {5}),
          ("P3", {5, 7}),
          ("P4", {2}),
          ("P5", {2, 7}),
          ("P6", {2, 5, 7}),
          ("P7", {1}),
          ("P8", {1, 5}),
          ("P9", {1, 2}),
          ("P10", {1, 2, 5, 7}),
          ("P11", {0}),
          ("P12", {0, 8}),
          ("P13", {0, 5}),
          ("P14", {0, 2}),
          ("P15", {0, 1, 5}),
          ("P16", {0, 1, 5, 7}),
          ("P17", {0, 1, 2, 5})],
        {0,1,2,5,7,8}
        ]
    ]

    def test_petric_function(self):
        for implicants, coverage in self.known_values:
            result = cora._find_irredundant_sums(implicants, coverage)
            self.assertEqual(76, len(result))
if __name__ == '__main__':
    unittest.main()
