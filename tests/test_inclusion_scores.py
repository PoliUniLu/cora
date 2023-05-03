import unittest
import pandas as pd
import cora
from cora import OptimizationContext
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def fpath(fname):
    return os.path.join(THIS_DIR, fname)


class TestInclusions(unittest.TestCase):
    def test_inclusion_score2(self):
        data = pd.read_csv(fpath("data/amenta_elliot_2017_5_causes_cs.csv"))
        context = OptimizationContext(data,["R","C"],case_col="CaseID",
                                      inc_score1=0.6,inc_score2=0.5,U=1)
        res1 = context.get_preprocessed_data()
        res = res1[["R", "C"]].to_dict(orient = 'list')
        expected_result = {'R':[1,1,1,1,0,1,0,0,1,1,0,1,1,1,1,0],
                           'C':[0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1]}

        self.assertEqual(res,expected_result)



