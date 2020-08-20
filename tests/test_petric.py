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
            
    def test_prepare_rows(self):
        data = pandas.read_json('{"ZAJ":{"0":0,"1":0,"2":0,"3":0,"4":1,"5":1,"6":1,"7":1,"8":1,"9":1,"10":1},"IC":{"0":0,"1":0,"2":1,"3":1,"4":0,"5":0,"6":0,"7":1,"8":1,"9":1,"10":1},"RES":{"0":2,"1":3,"2":1,"3":4,"4":1,"5":2,"6":4,"7":1,"8":2,"9":3,"10":4},"EURO":{"0":0,"1":0,"2":1,"3":1,"4":0,"5":0,"6":1,"7":0,"8":0,"9":1,"10":1},"E":{"0":1,"1":1,"2":1,"3":0,"4":1,"5":0,"6":1,"7":1,"8":0,"9":0,"10":1}}')
        outcol = ["E"]
        processed = cora.prepareRows(data, outcol)
        print(processed)
        self.assertTrue(processed is not None)

if __name__ == '__main__':
    unittest.main()
