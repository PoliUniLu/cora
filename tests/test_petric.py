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
            result=cora.find_irredundant_sums(implicants,coverage)
            self.assertEqual(irr_sums,result)

if __name__ == '__main__':
    unittest.main()
