"""Unit tests for benchmark.parser module."""

import pytest

from benchmark.parser import parse_jvalues
from benchmark.types import JValueBlock


class TestParseJvaluesBasic:
    """Tests for basic parse_jvalues functionality."""

    def test_empty_input(self):
        """Empty string returns no blocks."""
        result = parse_jvalues("")
        assert result == []

    def test_no_marker(self):
        """Input without marker returns no blocks."""
        result = parse_jvalues("Some random output\nNo J-values here\n")
        assert result == []

    def test_single_block_cpp_format(self):
        """Parse a single block in C++ output format (3 species)."""
        stdout = """\
some debug output
 Fast-J ----J-values----
L=      O2       O3   O3(1D)   
 57 0.00E+00 0.00E+00 0.00E+00
 56 2.54E-09 8.71E-03 7.38E-03
 55 1.34E-09 8.56E-03 7.24E-03
other output here
"""
        result = parse_jvalues(stdout)
        assert len(result) == 1
        block = result[0]
        assert block.sza_index == 0
        assert block.species == ["O2", "O3", "O3(1D)"]
        assert len(block.values) == 3
        assert block.values[57] == [0.0, 0.0, 0.0]
        assert block.values[56] == pytest.approx([2.54e-09, 8.71e-03, 7.38e-03])
        assert block.values[55] == pytest.approx([1.34e-09, 8.56e-03, 7.24e-03])

    def test_single_block_fortran_format(self):
        """Parse a single block in Fortran output format (many species)."""
        stdout = """\
   1   0.05   0.00   0.05
 Fast-J ----J-values----
 L=  O2       O3       O3(1D)   NO       H2COa
 57 2.54E-09 8.71E-03 7.38E-03 2.09E-06 1.03E-04
 56 1.34E-09 8.56E-03 7.24E-03 1.53E-06 1.02E-04
some trailing output
"""
        result = parse_jvalues(stdout)
        assert len(result) == 1
        block = result[0]
        assert block.sza_index == 0
        assert block.species == ["O2", "O3", "O3(1D)", "NO", "H2COa"]
        assert len(block.values) == 2
        assert block.values[57] == pytest.approx(
            [2.54e-09, 8.71e-03, 7.38e-03, 2.09e-06, 1.03e-04]
        )
        assert block.values[56] == pytest.approx(
            [1.34e-09, 8.56e-03, 7.24e-03, 1.53e-06, 1.02e-04]
        )

    def test_multiple_blocks(self):
        """Parse multiple SZA blocks from a single execution."""
        stdout = """\
debug output
 Fast-J ----J-values----
L=      O2       O3   O3(1D)   
 57 0.00E+00 0.00E+00 0.00E+00
  1 0.00E+00 0.00E+00 0.00E+00
other stuff between blocks
more debug
 Fast-J ----J-values----
L=      O2       O3   O3(1D)   
 57 1.00E-10 2.00E-10 3.00E-10
  1 4.00E-10 5.00E-10 6.00E-10
trailing output
 Fast-J ----J-values----
L=      O2       O3   O3(1D)   
 57 7.00E-10 8.00E-10 9.00E-10
  1 1.00E-09 1.10E-09 1.20E-09
end
"""
        result = parse_jvalues(stdout)
        assert len(result) == 3
        assert result[0].sza_index == 0
        assert result[1].sza_index == 1
        assert result[2].sza_index == 2
        # Check second block values
        assert result[1].values[57] == pytest.approx([1.0e-10, 2.0e-10, 3.0e-10])
        assert result[1].values[1] == pytest.approx([4.0e-10, 5.0e-10, 6.0e-10])
        # Check third block values
        assert result[2].values[57] == pytest.approx([7.0e-10, 8.0e-10, 9.0e-10])

    def test_variable_species_count(self):
        """Handle different number of species columns between blocks."""
        stdout = """\
 Fast-J ----J-values----
L=      O2       O3   O3(1D)   
  3 1.00E-05 2.00E-05 3.00E-05
  2 4.00E-05 5.00E-05 6.00E-05
  1 7.00E-05 8.00E-05 9.00E-05
"""
        result = parse_jvalues(stdout)
        assert len(result) == 1
        assert result[0].species == ["O2", "O3", "O3(1D)"]
        assert len(result[0].values) == 3
        assert result[0].values[3] == pytest.approx([1.0e-05, 2.0e-05, 3.0e-05])
        assert result[0].values[2] == pytest.approx([4.0e-05, 5.0e-05, 6.0e-05])
        assert result[0].values[1] == pytest.approx([7.0e-05, 8.0e-05, 9.0e-05])

    def test_single_digit_layer_numbers(self):
        """Layer numbers with varying widths are parsed correctly."""
        stdout = """\
 Fast-J ----J-values----
L=      O2       O3   
 10 1.00E-05 2.00E-05
  9 3.00E-05 4.00E-05
  1 5.00E-05 6.00E-05
"""
        result = parse_jvalues(stdout)
        assert len(result) == 1
        assert 10 in result[0].values
        assert 9 in result[0].values
        assert 1 in result[0].values


class TestParseJvaluesEdgeCases:
    """Edge case tests for parse_jvalues."""

    def test_marker_at_end_of_output(self):
        """Marker at end of file with no following content."""
        stdout = "some output\n Fast-J ----J-values----\n"
        result = parse_jvalues(stdout)
        assert result == []

    def test_marker_without_header(self):
        """Marker followed by non-header line is skipped."""
        stdout = """\
 Fast-J ----J-values----
This is not a header line
 57 1.00E-05 2.00E-05
"""
        result = parse_jvalues(stdout)
        # Should skip this block since no "L=" header found
        assert result == []

    def test_zero_values(self):
        """All-zero blocks are parsed correctly."""
        stdout = """\
 Fast-J ----J-values----
L=      O2       O3   O3(1D)   
 57 0.00E+00 0.00E+00 0.00E+00
 56 0.00E+00 0.00E+00 0.00E+00
"""
        result = parse_jvalues(stdout)
        assert len(result) == 1
        assert result[0].values[57] == [0.0, 0.0, 0.0]
        assert result[0].values[56] == [0.0, 0.0, 0.0]

    def test_large_and_small_exponents(self):
        """Values with various exponent magnitudes are parsed correctly."""
        stdout = """\
 Fast-J ----J-values----
L=      O2       O3   
  5 1.23E+02 4.56E-22
  4 7.89E+00 1.01E-01
"""
        result = parse_jvalues(stdout)
        assert len(result) == 1
        assert result[0].values[5] == pytest.approx([1.23e+02, 4.56e-22])
        assert result[0].values[4] == pytest.approx([7.89e+00, 1.01e-01])

    def test_block_terminated_by_empty_line(self):
        """A block terminated by an empty line is parsed correctly."""
        stdout = """\
 Fast-J ----J-values----
L=      O2       O3   
  2 1.00E-05 2.00E-05
  1 3.00E-05 4.00E-05

Some other output after empty line
"""
        result = parse_jvalues(stdout)
        assert len(result) == 1
        assert len(result[0].values) == 2


class TestParseJvaluesRealOutput:
    """Integration test using actual reference output format."""

    def test_parse_cpp_reference_format(self):
        """Verify parsing of actual C++ reference output structure."""
        # Simplified version of actual cpp_reference_output.txt
        stdout = """\
l2lev[57] = 57, jx=0, l2lev[l-1]=56
l2lev[58] = 58, jx=0, l2lev[l-1]=57
 Fast-J ----J-values----
L=      O2       O3   O3(1D)   
 57 0.00E+00 0.00E+00 0.00E+00
 56 0.00E+00 0.00E+00 0.00E+00
 55 0.00E+00 0.00E+00 0.00E+00
 41 4.78E-21 4.78E-20 8.53E-21
 40 2.02E-20 2.02E-19 3.03E-20
  2 0.00E+00 0.00E+00 0.00E+00
  1 0.00E+00 0.00E+00 0.00E+00
l2lev[1] = 1, jx=0, l2lev[l-1]=0
l2lev[2] = 2, jx=0, l2lev[l-1]=1
"""
        result = parse_jvalues(stdout)
        assert len(result) == 1
        block = result[0]
        assert block.species == ["O2", "O3", "O3(1D)"]
        assert block.values[41] == pytest.approx([4.78e-21, 4.78e-20, 8.53e-21])
        assert block.values[1] == [0.0, 0.0, 0.0]
        # Block ends at non-numeric line "l2lev..."
        assert len(block.values) == 7
