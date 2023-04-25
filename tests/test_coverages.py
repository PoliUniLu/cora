import unittest
import pandas as pd
import os
from cora import OptimizationContext,data_mining
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

def fpath(fname):
    return os.path.join(THIS_DIR, fname)
class OptimizationContext_tests(unittest.TestCase):


    def test_coverage_scores_multi_output(self):
        data = pd.read_csv(fpath('data/ansorg_2014.csv'))
        context = OptimizationContext(data,["INTGOV","ENGIC"],
                                      input_labels=["ECOFAIL",
                                                    "STATEFAIL",
                                                    "REGID",
                                                    "ECONET"],inc_score1=0.5)
        result_pi_details = context.pi_details()
        res_final_pi_details = [v.dropna().to_dict()
                                for k, v in result_pi_details.iterrows()]
        expected_result_pi_details = [{'PI': 'ECOFAIL{1}*ECONET{1}',
                            'Cov.r': 0.2, 'Inc.': 0.5, 'S1': 0.2, 'S5': 0.2,
                            'S6': 0.2, 'S8': 0.2, 'S9': 0.2, 'S14': 0.2},
                           {'PI': 'STATEFAIL{1}*ECONET{1}',
                            'Cov.r': 0.2, 'Inc.': 0.5, 'S3': 0.2, 'S4': 0.2,
                            'S7': 0.2, 'S11': 0.2, 'S12': 0.2, 'S15': 0.2},
                           {'PI': 'REGID{1}*ECONET{1}', 'Cov.r': 0.4,
                            'Inc.': 0.67, 'S2': 0.4, 'S10': 0.4, 'S13': 0.4},
                           {'PI': 'REGID{1}*ECONET{2}', 'Cov.r': 0.6,
                            'Inc.': 0.6, 'S1': 0.6, 'S2': 0.6, 'S3': 0.6,
                            'S4': 0.6, 'S5': 0.6, 'S6': 0.6, 'S7': 0.6, ''
                            'S8': 0.6, 'S9': 0.6, 'S10': 0.6, 'S11': 0.6,
                            'S12': 0.6, 'S13': 0.6, 'S14': 0.6, 'S15': 0.6},
                           {'PI': 'ECOFAIL{1}*REGID{0}', 'Cov.r': 0.5,
                            'Inc.': 1.0, 'S2': 0.5, 'S3': 0.5, 'S9': 0.5,
                            'S12': 0.5, 'S14': 0.5},
                           {'PI': 'STATEFAIL{0}*REGID{1}', 'Cov.r': 0.2,
                            'Inc.': 1.0, 'S1': 0.2, 'S3': 0.2, 'S7': 0.2,
                            'S8': 0.2, 'S9': 0.2, 'S11': 0.2},
                           {'PI': 'STATEFAIL{2}*REGID{1}', 'Cov.r': 0.2,
                            'Inc.': 0.5},
                           {'PI': 'ECOFAIL{0}*REGID{1}',
                                           'Cov.r': 0.2, 'Inc.': 1.0,
                                           'S4': 0.2, 'S5': 0.2, 'S6': 0.2,
                            'S12': 0.2, 'S14': 0.2, 'S15': 0.2},
                           {'PI': 'ECONET{2}', 'Cov.r': 1.0, 'Inc.': 1.0,
                            'S1': 0.5, 'S4': 0.5, 'S6': 0.5,
                            'S10': 0.5, 'S11': 0.5},
                           {'PI': 'STATEFAIL{2}', 'Cov.r': 0.5,
                            'Inc.': 1.0, 'S5': 0.25, 'S7': 0.25, 'S8': 0.25,
                            'S13': 0.25, 'S15': 0.25},
                           {'PI': 'STATEFAIL{1}*REGID{0}', 'Cov.r': 0.25,
                            'Inc.': 1.0, 'S5': 0.25, 'S7': 0.25, 'S8': 0.25,
                            'S13': 0.25, 'S15': 0.25}]
        result_system_details = context.system_details()
        expected_result_system_details = {'Cov.': {'Solution details': 1.0},
                                          'Inc.': {'Solution details': 0.9}}
        self.assertEqual(result_system_details.to_dict(),
                         expected_result_system_details)
        #self.assertEqual(res_final_pi_details,expected_result_pi_details)

    def test_coverage_scores_single_output(self):
         df = pd.DataFrame([[1,1,0,1],
                            [0,0,1,1],
                            [1,0,1,0],
                            [0,1,0,1]],
                        columns=["A","B","C","OUT"])
         context = OptimizationContext(data = df,
                                      output_labels = ["OUT"])
         result = context.pi_details()
         res_final = [v.dropna().to_dict() for k, v in result.iterrows()]
         expected_result = [{'PI': 'B', 'Cov.r': 0.67, 'Inc.': 1.0, 'M1': 0.33},
                            {'PI': 'c', 'Cov.r': 0.67, 'Inc.': 1.0, 'M2': 0.33},
                            {'PI': '#a', 'Cov.r': 0.67, 'Inc.': 1.0, 'M1': 0.33,
                             'M2': 0.33}]
         self.assertEqual(res_final,expected_result)
    def test_system_details(self):
        df = pd.DataFrame([[1, 1, 0, 1],
                           [0, 0, 1, 1],
                           [1, 0, 1, 0],
                           [0, 1, 0, 1]],
                          columns = ["A", "B", "C", "OUT"])
        context = OptimizationContext(data=df,
                                            output_labels = ["OUT"])
        result = context.system_details()
        expected_result = {'Cov.': {'Solution details': 1.0},
                           'Inc.': {'Solution details': 1.0}}
        self.assertEqual(result.to_dict(),expected_result)



