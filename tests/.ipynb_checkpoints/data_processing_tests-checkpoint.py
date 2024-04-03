# +
from qPCR_analysis import Data_processing
import pytest 

def test_primer_efficiency_calc():
    assert Data_processing.primer_effiency_calc(-3.754) == pytest.approx(84.66, 0.01), "Primer efficiency calculation with a slope of -3.754 should equal 84.66 " 

def test_delta_Ct():
    assert Data_processing.delta_Ct(19,18) == pytest.approx(1), "delta_Ct calculation with a Ct of 19 and 18 should be 1" 
    
def test_pfaffl_calc():
    assert Data_processing.pfaffl_calc(100, 2, 5, 2) == pytest.approx(400), "pfaffl_calc with efficiency of 100 and 5 and delta ct of 2 and 2 should be 400"

def test_calculate_polysome_diff(reference_cts, cts):
    assert Data_processing.delta_Ct(19,18) == pytest.approx(1), "test_calculate_polysome_diff calculation with a Ct of 19 and 18 should be 1" 
    
def test_calculate_2_delta_ct_polysome(delta_ct):
    assert Data_processing.calculate_2_delta_ct_polysome(2) == pytest.approx(4)

def test_calculate_percentages(total, column):
    df = pd.Data
    assert Data_processing.calculate_percentages(
    
    return (column * 100) / total if total != 0 else 0
