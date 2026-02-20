import unittest
import pandas as pd
from parameterized import parameterized
from cora import OptimizationContext
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def fpath(fname):
    return os.path.join(THIS_DIR, fname)


class TestPrimeImplicants(unittest.TestCase):

    INPUT2 = pd.read_csv(fpath("data/Multi_output/data2.csv"), sep=",")
    INPUT3 = pd.read_csv(fpath("data/Multi_output/data3.csv"), sep=",")
    INPUT5 = pd.read_csv(fpath("data/Multi_output/data5.csv"), sep=",")
    INPUT6 = pd.read_csv(fpath("data/Multi_output/data6.csv"), sep=",")

    @parameterized.expand(
        [
            (
                INPUT2,
                (
                    "x*y*F1",
                    "y*z*F1",
                    "W*X*F1",
                    "W*y*F1",
                    "W*X*Z",
                    "W*x*y",
                    "w*X*Z",
                    "W*y*Z",
                    "w*F1",
                    "X*Y",
                    "Y*f1",
                    "Y",
                ),
            ),
            (
                INPUT3,
                (
                    (
                        "a*b*c*f1",
                        "A*B*C*d",
                        "C*d*F1",
                        "a*b*c",
                        "a*c*F1",
                        "a*c*D",
                        "b*d*F1",
                        "b*c*d",
                        "A*b*F1",
                        "C*d",
                        "b*d",
                        "A*f1",
                        "b*F1",
                        "D*f1",
                        "b*f1",
                        "C*f1",
                        "A*b",
                        "b*C",
                    )
                ),
            ),
            (
                INPUT5,
                (
                    "a*C*d*f1",
                    "a*B*C",
                    "c*d*F1",
                    "a*c*F1",
                    "a*D*F1",
                    "A*b*d",
                    "A*C*F1",
                    "b*C*d",
                    "a*b*F1",
                    "A*c*d",
                    "A*b*C",
                    "B*c*F1",
                    "b*C*D",
                    "B*c*d",
                    "a*C*D",
                    "A*b*f1",
                    "b*C*f1",
                    "C*F1",
                ),
            ),
            (
                INPUT6,
                (
                    "A*b*d",
                    "a*C*D",
                    "a*B*D",
                    "b*C*d",
                    "a*b*C",
                    "a*D",
                    "a*C",
                    "a*b",
                    "b*d",
                ),
            ),
        ]
    )
    def test_prime_implicants_multiple(self, data, expected_pis):
        prime_impl = OptimizationContext(
            data, output_labels=list(data.columns[-2:])
        ).get_prime_implicants()
        res = tuple(str(i) for i in prime_impl)
        self.assertEqual(res, expected_pis)

    INPUT1 = pd.read_csv(fpath("data/One_output/data1.csv"), sep=",")
    INPUT2 = pd.read_csv(fpath("data/One_output/data2.csv"), sep=",")
    INPUT3 = pd.read_csv(fpath("data/One_output/data3.csv"), sep=",")
    INPUT4 = pd.read_csv(fpath("data/One_output/data4.csv"), sep=",")
    INPUT5 = pd.read_csv(fpath("data/One_output/data5.csv"), sep=",")

    @parameterized.expand(
        [
            (INPUT1, {"a*c*D", "a*B*C", "a*B*D", "#b*c", "#C*d"}),
            (INPUT2, {"a*b", "#C*D", "#A*D", "b*C"}),
            (
                INPUT3,
                {
                    "#v*W*x*Y",
                    "#W*X*Y*z",
                    "#W*X*y*Z",
                    "#V*w*x*Y",
                    "#v*w*X*Y",
                    "#w*x*z",
                    "#v*z",
                },
            ),
            (INPUT4, {"#C*D", "#A*D", "#b*C"}),
            (
                INPUT5,
                {
                    "v*w*X*Z",
                    "v*w*Y*Z",
                    "v*W*y*Z",
                    "v*w*x*Y",
                    "v*X*y*Z",
                    "v*W*x*Z",
                    "v*x*Y*Z",
                    "#V*W*z",
                    "#w*x*z",
                    "#W*X*Y*z",
                },
            ),
        ]
    )
    def test_prime_implicants_single(self, data, expected_pis):
        prime_impl_on_off = OptimizationContext(
            data, output_labels=list(data.columns[-1]), algorithm="ON-OFF"
        ).get_prime_implicants()
        prime_impl_on_dc = OptimizationContext(
            data,
            output_labels=list(data.columns[-1]),
        ).get_prime_implicants()
        res_on_off = set(str(i) for i in prime_impl_on_off)
        res_on_dc = set(str(i) for i in prime_impl_on_dc)
        self.assertEqual(res_on_dc, res_on_off)
        self.assertEqual(expected_pis, res_on_off)
        self.assertEqual(expected_pis, res_on_dc)


if __name__ == "__main__":
    unittest.main()
