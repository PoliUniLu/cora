import cora
import pandas as pd
from cora import OptimizationContext
import unittest
from pandas.testing import assert_frame_equal


class TestOutputs(unittest.TestCase):
    def test_output_exeption(self):
        df = pd.DataFrame(
            [[1, 0, 1, 0], [0, 1, 1, 1], [0, 1, 0, 0], [1, 0, 1, 2]],
            columns=["a", "b", "c", "o"],
        )
        context = cora.OptimizationContext(df, ["o"])
        with self.assertRaises(Exception):
            preprocessed_data = context._preprocess_data()

    def test_output(self):
        df = pd.DataFrame(
            [[1, 0, 1, 1], [0, 1, 1, 1], [0, 1, 0, 0], [1, 0, 1, 2]],
            columns=["a", "b", "c", "o"],
        )
        context1 = cora.OptimizationContext(df, ["o{1}"])
        df_test1 = pd.DataFrame(
            [[0, 1, 0, 0], [0, 1, 1, 1], [1, 0, 1, 0]], columns=["a", "b", "c", "o"]
        )
        df_test1.astype({"o": "int64"})
        assert_frame_equal(df_test1, context1.get_preprocessed_data())
        context2 = cora.OptimizationContext(df, ["o{1,2}"])
        df_test2 = pd.DataFrame(
            [[0, 1, 0, 0], [0, 1, 1, 1], [1, 0, 1, 1]], columns=["a", "b", "c", "o"]
        )
        df_test2.astype({"o": "int64"})
        assert_frame_equal(df_test2, context2.get_preprocessed_data())


if __name__ == "__main__":
    unittest.main()
