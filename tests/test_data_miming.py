import unittest
import pandas as pd
from cora import data_mining
from pandas.testing import assert_frame_equal
import os
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

def fpath(fname):
    return os.path.join(THIS_DIR, fname)


class TestPrimeImplicants(unittest.TestCase):

    def test_data_mining(self):
        data = pd.DataFrame([[1,0,1,0],[1,0,0,0],[0,1,0,1],[0,0,1,1]],
                            columns = ["A","B","C","O"])
        data_mining_on_off = data_mining(data,["O"],3,algorithm='ON-OFF')
        data_mining_on_dc = data_mining(data,["O"],3,algorithm='ON-DC')

        assert_frame_equal(data_mining_on_dc,data_mining_on_off)

    def test_data_mining_multi_output(self):
        data = pd.DataFrame([[1,0,1,0,1],[1,0,0,0,1],[0,1,0,1,1],[0,0,1,1,0]],
                            columns = ["A","B","C","Y","X"])

        data_mining_on_off = data_mining(data,["Y","X"],1,algorithm="ON-OFF")
        print(data_mining_on_off)
        data_mining_on_dc = data_mining(data, ["Y","X"], 1,algorithm="ON-DC")

        assert_frame_equal(data_mining_on_dc,data_mining_on_off)
    def test_data_mining_single_output(self):
        data =  pd.read_csv(fpath('data/data_mining_tests.csv'))
        print(data)

        res_on_off = data_mining(data, ["XWES3{1}"], 2,
                          input_labels=["XNS3","XPS3","XGHG3"],
                          case_col='COUNTRY',
                          inc_score1=0.5,
                          algorithm="ON-OFF")
        res_on_dc = data_mining(data, ["XWES3{1}"], 2,
                          input_labels=["XNS3","XPS3","XGHG3"],
                          case_col='COUNTRY',
                          inc_score1=0.5,
                          algorithm="ON-DC")
        assert_frame_equal(res_on_dc, res_on_off)
    def test_data_mining_automatic(self):
        data = pd.read_csv(fpath('data/testdata1.csv'))
        res_on_off_automatic = data_mining(data, ["S"], 1,
                                 input_labels=["C",  "R" ,"L"  ,"E"],
                                 case_col='Case',
                                 inc_score1=0.8,
                                 algorithm="ON-OFF",
                                 automatic = True)
        res_on_dc_nonautomatic = data_mining(data, ["S"], 2,
                                 input_labels=["C",  "R" ,"L"  ,"E"],
                                 case_col='Case',
                                 inc_score1=0.8,
                                 algorithm="ON-DC"
                                 )

        assert_frame_equal(res_on_dc_nonautomatic, res_on_off_automatic)



if __name__ == '__main__':
    unittest.main()
