import unittest
import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from cora import OptimizationContextV2


class PreprocessingTest(unittest.TestCase):

    def test_preprocess_empty(self):
        opt_context = OptimizationContextV2(pd.DataFrame(), [])
        result = opt_context.preprocess_data()
        self.assertTrue(result.empty)

    # test duplicities
    def test_preprocess_simple(self):
        opt_context = OptimizationContextV2(pd.DataFrame(data=[[1, 1, 1, 1],
                                                               [1, 1, 1, 1],
                                                               [1, 0, 0, 0]],
                                                         columns=["A", "B", "C", "O"]),
                                            ["O"])
        result = opt_context.preprocess_data()

        self.assertEqual(len(result), 2)

    def test_inclusion_scores(self):
        opt_context = OptimizationContextV2(pd.DataFrame(data=[["A", 1, 1, 1, 0],
                                                               ["B", 1, 1, 1, 1],
                                                               ["C", 1, 0, 0, 0]],
                                                         columns=["ID", "A", "B", "C", "O"]),
                                            ["O"],
                                            case_col="ID",
                                            inc_score1=0.6)
        result = opt_context.preprocess_data()
        expected_cols = {
            "A": [1, 1],
            "B": [0, 1],
            "C": [0, 1],
            "n": [1, 2],
            "Cases": ["C", "A,B"],
            "Inc_O": [0.0, 0.5],
            "O": [0, 0]
        }
        expected_result = pd.DataFrame.from_dict(expected_cols)
        pd.testing.assert_frame_equal(result, expected_result)

    def test_levels(self):
        opt_context = OptimizationContextV2(pd.DataFrame(data=[[ 1, 1, 1, 0],
                                                               [ 1, 1, 1, 1],
                                                               [ 1, 0, 0, 0]],
                                                         columns=["A", "B", "C", "O"]),
                                           ["O"])
        result = opt_context.get_levels()
        print(result)
        self.assertEqual(result,([2,2,2],[[0,1],[0,1],[0,1]]))

    def test_prepare_rows(self):
        opt_context = OptimizationContextV2(pd.DataFrame(data=[[1, 1, 1, 0],
                                                               [1, 1, 1, 1],
                                                               [1, 0, 0, 0]],
                                                         columns=["A", "B", "C", "O"]),
                                            ["O"])
        result = opt_context.prepareRows()
        expected_result = np.array([[0,0,0],
                            [0,0,1],
                            [0,1,0],
                            [0,1,1],
                            [1,0,1],
                            [1,1,0]])

        np.testing.assert_array_equal(result, expected_result)




if __name__ == '__main__':
    unittest.main()
