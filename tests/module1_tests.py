# +
from qPCR_analysis import module1
import pytest 

def test_primer_efficiency_calc():
    # Assuming the function returns a float and using `pytest.approx` for comparison with a tolerance
    assert module1.primer_effiency_calc(-3.754) == pytest.approx(84.66, 0.01)


