import unittest
import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0,
                os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from cora import OptimizationContext

# test the Class Optimization Context
class OptimizationContext_tests(unittest.TestCase):

    def test_inclusion_scores(self):
        opt_context = OptimizationContext(pd.DataFrame(data=[["1", 1, 0, 1, 0],
                                                             ["2", 1, 1, 1, 0],
                                                             ["3", 0, 1, 0, 1],
                                                             ["4", 1 ,1, 1, 1],
                                                             ["5", 1, 0, 1, 1]],
                                                       columns=["ID", "A", "B",
                                                                "C","O"]),
                                          ["O"],
                                          case_col="ID",
                                          inc_score1=0.5)
        opt_context._preprocess_data()
        result = opt_context.preprocessed_data_raw
        expected_cols = {
            "A": [0, 1, 1],
            "B": [1, 0, 1],
            "C": [0, 1, 1],
            "n": [1, 2, 2],
            "Cases": ["3", "1,5", "2,4"],
            "Inc_O": [1.0, 0.5, 0.5],
            "O": [1, 1, 1]
        }
        expected_result = pd.DataFrame.from_dict(expected_cols)
        pd.testing.assert_frame_equal(result, expected_result)

        opt_context = OptimizationContext(
                pd.DataFrame(data=[["1", 1, 0, 1, 0, 1],
                                   ["2", 1, 1, 1, 0, 1],
                                   ["3", 0, 1, 0, 1, 0],
                                   ["4", 1, 1, 1, 1, 1],
                                   ["5", 1, 0, 1, 1, 1]],
                             columns=["ID", "A", "B",
                                      "C", "O","P"]),
                ["O","P"],
                case_col="ID",
                inc_score1=0.5)
        opt_context._preprocess_data()
        result = opt_context.preprocessed_data_raw
        expected_cols = {
                "A": [0, 1, 1],
                "B": [1, 0, 1],
                "C": [0, 1, 1],
                "n": [1, 2, 2],
                "Cases": ["3", "1,5", "2,4"],
                "Inc_O": [1.0, 0.5, 0.5],
                "Inc_P": [0, 1, 1],
                "O": [1, 1, 1],
                "P": [0, 1, 1]

            }
        expected_result = pd.DataFrame.from_dict(expected_cols)
        pd.testing.assert_frame_equal(result, expected_result)

        opt_context = OptimizationContext(pd.DataFrame(data=[["1", 1, 0, 1, 0],
                                                             ["2", 1, 1, 1, 0],
                                                             ["3", 0, 1, 0, 1],
                                                             ["4", 1, 1, 1, 1],
                                                             ["5", 1, 0, 1, 1]],
                                                       columns=["ID", "A", "B",
                                                                "C", "O"]),
                                          ["O"],
                                          case_col="ID",
                                          inc_score1=0.6)
        opt_context._preprocess_data()
        result = opt_context.preprocessed_data_raw
        expected_cols = {
            "A": [0, 1, 1],
            "B": [1, 0, 1],
            "C": [0, 1, 1],
            "n": [1, 2, 2],
            "Cases": ["3", "1,5", "2,4"],
            "Inc_O": [1.0, 0.5, 0.5],
            "O": [1, 0, 0]
        }
        expected_result = pd.DataFrame.from_dict(expected_cols)
        pd.testing.assert_frame_equal(result, expected_result)


    def test_levels(self):
        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 1, 1, 1],
                                                             [0, 1, 1, 1],
                                                             [1, 0, 0, 1]],
                                                       columns=["A", "B", "C",
                                                                "O"]),
                                          ["O"])
        opt_context._preprocess_data()
        result = opt_context.get_levels()
        self.assertEqual(result, ([2, 2, 2]))

        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 3, 6, 1, 1],
                                                             [0, 1, 1, 1, 0],
                                                             [1, 3, 0, 1, 0]],
                                                       columns=["A", "B", "C",
                                                                "O","P"]),
                                          ["O","P"])
        opt_context._preprocess_data()
        result = opt_context.get_levels()
        self.assertEqual(result, ([2, 2, 3]))


        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 1, 1, 1],
                                                             [3, 1, 1, 1],
                                                             [2, 0, 2, 1]],
                                                       columns=["A", "B", "C",
                                                                "O"]),
                                          ["O"])
        opt_context._preprocess_data()
        result = opt_context.get_levels()
        self.assertEqual(result, ([3, 2, 2]))

    def test_prepare_rows(self):
        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 1, 1, 0],
                                                             [1, 1, 1, 1],
                                                             [0, 0, 0, 1]],
                                                       columns=["A", "B", "C",
                                                                "O"]),
                                          ["O"])
        opt_context._prepareRows()
        result = opt_context.table
        expected_result = np.array([[0, 0, 0],
                                    [0, 0, 1],
                                    [0, 1, 0],
                                    [0, 1, 1],
                                    [1, 0, 0],
                                    [1, 0, 1],
                                    [1, 1, 0]])


        np.testing.assert_array_equal(result, expected_result)


        np.testing.assert_array_equal(result, expected_result)
    def test_get_prime_implicants_1_DC(self):
        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 1, 1, 1],
                                                             [0, 1, 1, 1],
                                                             [1, 0, 0, 0]],
                                                       columns=["A", "B", "C",
                                                                "O"]),
                                          ["O"])
        result = opt_context.get_prime_implicants_1_DC()

        self.assertEqual(len(result), 3)
        self.assertEqual(set(str(impl) for impl in result), {'a', 'C', 'B'})

        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 2, 1, 1, 1],
                                                             [1, 2, 0, 1, 1],
                                                             [2, 1, 0, 0, 0]],
                                                       columns=["A", "B", "C",
                                                                "D",
                                                                "O"]),
                                          ["O"])
        result = opt_context.get_prime_implicants_1_DC()

        self.assertEqual(len(result), 4)
        self.assertEqual(set(str(impl) for impl in result),
                         {'C', 'A', 'D', 'B'})

        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 3, 1, 1, 1],
                                                             [1, 2, 0, 2, 1],
                                                             [2, 1, 0, 0, 0]],
                                                       columns=["A", "B", "C",
                                                                "D",
                                                                "O"]),
                                          ["O"])
        result = opt_context.get_prime_implicants_1_DC()

        self.assertEqual(len(result), 6)
        self.assertEqual(set(str(impl) for impl in result),
                    {'B{3}', 'D{2}', 'C{1}', 'D{1}', 'A{1}', 'B{2}'})

        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 2, 1, 1, 1],
                                                             [1, 1, 2, 1, 1],
                                                             [1, 2, 1, 2, 1],
                                                             [0, 1, 0, 0, 0],
                                                             [1, 1, 0, 1, 0],
                                                             [1, 1, 1, 1 ,0],
                                                             [0, 0, 0, 1, 0],
                                                             [1, 2, 0, 2, 1],
                                                             [2, 1, 0, 0, 0]],
                                                       columns=["A", "B", "C",
                                                                "D",
                                                                "O"]),
                                          ["O"])
        result = opt_context.get_prime_implicants_1_DC()

        self.assertEqual(len(result), 2)
        self.assertEqual(set(str(impl) for impl in result),
                         {'#C{2}', '#B{2}'})

        opt_context = OptimizationContext(pd.DataFrame(data=[[1, 2, 1, 1, 1, 1],
                                                             [1, 1, 2, 1, 1, 1],
                                                             [1, 2, 1, 2, 1, 0],
                                                             [0, 1, 0, 0, 0, 0],
                                                             [1, 1, 0, 1, 0, 1],
                                                             [1, 1, 1, 1, 0, 1],
                                                             [0, 0, 0, 1, 0, 1],
                                                             [1, 2, 0, 2, 1, 0],
                                                            [2, 1, 0, 0, 0, 1]],
                                                       columns=["A", "B", "C",
                                                                "D",
                                                                "O","P"]),
                                          ["O","P"])
        result = opt_context.get_prime_implicants_1_DC()

        self.assertEqual(len(result), 9)
        self.assertEqual(set(str(impl) for impl in result),
                         {'A{1}*B{1}','A{2}','B{0}','B{1}*C{1}','B{2}',
                          'B{2}*D{1}','C{2}','D{1}','D{2}'})

        opt_context = OptimizationContext(pd.DataFrame(data=[[0, 2, 1, 0, 1, 1],
                                                             [1, 1, 2, 1, 1, 0],
                                                             ],
                                                       columns=["A", "B", "C",
                                                                "D",
                                                                "O", "P"]),
                                          ["O", "P"])
        result = opt_context.get_prime_implicants_1_DC()

        self.assertEqual(len(result), 5)
        self.assertEqual(set(str(impl) for impl in result),
                         {'B', '1', 'C', 'a', 'd'})

        def test_pi_chart(self):
            opt_context = OptimizationContext(
                pd.DataFrame(data=[[ 1, 0, 1, 0],
                                   [ 1, 1, 1, 0],
                                   [ 0, 1, 0, 1],
                                   [ 1, 1, 1, 1],
                                   [ 1, 0, 1, 1]],
                             columns=["A", "B",
                                      "C", "O"]),
                ["O"])

            result = opt_context.get_pi_chart()
            expected_cols = {
                "A": [0, 1, 1],
                "B": [1, 0, 1],
                "C": [0, 1, 1],
                "n": [1, 2, 2]

            }
            expected_result = pd.DataFrame.from_dict(expected_cols)
            pd.testing.assert_frame_equal(result, expected_result)


if __name__ == '__main__':
    unittest.main()
