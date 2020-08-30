import pytest

def test_convert_pvalue_to_text():

    p = 0.0
    bootstrap_replicates = 100000
    result = convert_pvalue_to_text(p, bootstrap_replicates)
    assert result == "<0.00001"

    p = 0.123
    bootstrap_replicates = 100000
    result = convert_pvalue_to_text(p, bootstrap_replicates)
    assert result == "0.12300"

    p = 0.001
    bootstrap_replicates = 100000
    result = convert_pvalue_to_text(p, bootstrap_replicates)
    assert result == "0.00100"

    p = 0.01234567890123456789
    bootstrap_replicates = 100000
    result = convert_pvalue_to_text(p, bootstrap_replicates)
    assert result == "0.01234"
