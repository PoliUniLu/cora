import unittest
import pandas as pd
from parameterized import parameterized
from cora import experimental_find_implicants_cubes



class TestPrimeImplicants(unittest.TestCase):

    INPUT1 = pd.DataFrame([[1,0,1,1],
                           [0,1,0,0],
                           [0,0,1,1],
                           [1,1,1,1]],
                          columns=["A","B","C","O"])
    INPUT2 = pd.DataFrame([[0,0,0,1],
                           [1,0,1,0],
                           [1,1,0,1],
                           [1,1,1,0]],
                          columns=["A","B","C","O"])

    @parameterized.expand([(INPUT1, [[1, -1, -1], [-1, 0, -1], [-1, -1, 1]]),
                           (INPUT2, [[0, -1, -1], [-1, -1, 0]])])
    def test_prime_implicants_single(self, data, expected_pis):

        prime_impl_cubes = experimental_find_implicants_cubes(
            data,
            outputs=list(data.columns[-1]),
            inputs = list(data.columns[:-1]))
        self.assertEqual(expected_pis, prime_impl_cubes)


if __name__ == '__main__':
    unittest.main()
