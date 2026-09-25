"""Property-based tests for benchmark.parser module.

**Validates: Requirements 4.1**

Property 2: J-value parsing round-trip consistency
For any well-formed J-value output block, `parse_jvalues` extracts all
layer/species values correctly.
"""

from hypothesis import given, settings, assume
from hypothesis import strategies as st

from benchmark.parser import parse_jvalues


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

# Strategy for species names: short alphanumeric identifiers, may include
# parentheses and digits like "O3(1D)", "H2COa"
species_name_st = st.from_regex(r"[A-Z][A-Za-z0-9()]{0,7}", fullmatch=True)

# Strategy for layer numbers (positive integers, typically 1-57)
layer_number_st = st.integers(min_value=1, max_value=99)

# Strategy for J-values: non-negative floats in scientific notation range
# Fortran e9.2 format covers roughly 0.00E+00 to 9.99E+99
jvalue_st = st.one_of(
    st.just(0.0),
    st.floats(min_value=1e-30, max_value=9.99e+30, allow_nan=False, allow_infinity=False),
)


def format_jvalue(val: float) -> str:
    """Format a float value in Fortran-like scientific notation (e.g., 2.54E-09)."""
    if val == 0.0:
        return "0.00E+00"
    # Use Python's %E format which is compatible with Fortran E notation
    return f"{val:.2E}"


def build_jvalue_block(species: list, layers: list, values: list) -> str:
    """Build a well-formed J-value output block string.

    Args:
        species: List of species name strings.
        layers: List of integer layer numbers (must be unique).
        values: 2D list [layer_idx][species_idx] of float values.

    Returns:
        A string representing a complete J-value output block.
    """
    lines = []
    lines.append(" Fast-J ----J-values----")

    # Header line: "L=" followed by species names
    header = "L=" + "".join(f"{s:>9s}" for s in species)
    lines.append(header)

    # Numeric rows: layer_num followed by values in scientific notation
    for layer_idx, layer_num in enumerate(layers):
        row_vals = " ".join(format_jvalue(values[layer_idx][s_idx]) for s_idx in range(len(species)))
        lines.append(f"{layer_num:3d} {row_vals}")

    # Empty line terminates the block
    lines.append("")
    return "\n".join(lines)


# Strategy for a complete well-formed J-value block configuration
@st.composite
def jvalue_block_config_st(draw):
    """Generate a valid J-value block configuration.

    Returns a tuple of (species, layers, values) that can be used to
    build a well-formed block and verify parsing results.
    """
    num_species = draw(st.integers(min_value=1, max_value=10))
    num_layers = draw(st.integers(min_value=1, max_value=20))

    species = draw(
        st.lists(
            species_name_st,
            min_size=num_species,
            max_size=num_species,
            unique=True,
        )
    )

    layers = draw(
        st.lists(
            layer_number_st,
            min_size=num_layers,
            max_size=num_layers,
            unique=True,
        )
    )

    # Generate values: one list per layer, one value per species
    values = []
    for _ in range(num_layers):
        row = draw(
            st.lists(jvalue_st, min_size=num_species, max_size=num_species)
        )
        values.append(row)

    return species, layers, values


# Strategy for multiple blocks (simulating multiple SZA outputs)
@st.composite
def multi_block_config_st(draw):
    """Generate 1-3 J-value block configurations."""
    num_blocks = draw(st.integers(min_value=1, max_value=3))
    configs = []
    for _ in range(num_blocks):
        config = draw(jvalue_block_config_st())
        configs.append(config)
    return configs


# ---------------------------------------------------------------------------
# Property 2: J-value parsing round-trip consistency
# ---------------------------------------------------------------------------


class TestJValueParsingProperty:
    """**Validates: Requirements 4.1**"""

    @given(config=jvalue_block_config_st())
    @settings(max_examples=200)
    def test_all_layer_numbers_extracted(self, config):
        """All layer numbers from the input are present in parsed output."""
        species, layers, values = config

        stdout = build_jvalue_block(species, layers, values)
        result = parse_jvalues(stdout)

        assert len(result) == 1
        block = result[0]
        assert set(block.values.keys()) == set(layers)

    @given(config=jvalue_block_config_st())
    @settings(max_examples=200)
    def test_all_species_names_extracted(self, config):
        """All species names from the header are extracted correctly."""
        species, layers, values = config

        stdout = build_jvalue_block(species, layers, values)
        result = parse_jvalues(stdout)

        assert len(result) == 1
        block = result[0]
        assert block.species == species

    @given(config=jvalue_block_config_st())
    @settings(max_examples=200)
    def test_values_parsed_within_precision(self, config):
        """All numeric values are parsed to within floating-point precision.

        Since we format with %.2E (3 significant digits) and parse back,
        the relative error should be within the representation precision
        of the e9.2 format.
        """
        species, layers, values = config

        stdout = build_jvalue_block(species, layers, values)
        result = parse_jvalues(stdout)

        assert len(result) == 1
        block = result[0]

        for layer_idx, layer_num in enumerate(layers):
            parsed_row = block.values[layer_num]
            assert len(parsed_row) == len(species)

            for s_idx in range(len(species)):
                original = values[layer_idx][s_idx]
                parsed = parsed_row[s_idx]

                if original == 0.0:
                    assert parsed == 0.0
                else:
                    # The e9.2 format gives ~3 significant digits, so
                    # relative error should be within 0.5%
                    rel_error = abs(parsed - original) / abs(original)
                    assert rel_error < 0.006, (
                        f"Layer {layer_num}, species {species[s_idx]}: "
                        f"original={original}, parsed={parsed}, rel_error={rel_error}"
                    )

    @given(config=jvalue_block_config_st())
    @settings(max_examples=200)
    def test_correct_number_of_values_per_layer(self, config):
        """Each layer row has exactly as many values as there are species."""
        species, layers, values = config

        stdout = build_jvalue_block(species, layers, values)
        result = parse_jvalues(stdout)

        assert len(result) == 1
        block = result[0]

        for layer_num in layers:
            assert len(block.values[layer_num]) == len(species)

    @given(configs=multi_block_config_st())
    @settings(max_examples=100)
    def test_multiple_blocks_all_parsed(self, configs):
        """Multiple J-value blocks are all parsed with correct sza_index."""
        # Build full stdout with multiple blocks and interleaving text
        parts = ["some preamble output\n"]
        for config in configs:
            species, layers, values = config
            parts.append(build_jvalue_block(species, layers, values))
            parts.append("some interleaving text\n")

        stdout = "\n".join(parts)
        result = parse_jvalues(stdout)

        assert len(result) == len(configs)

        for idx, (block, config) in enumerate(zip(result, configs)):
            species, layers, values = config
            assert block.sza_index == idx
            assert block.species == species
            assert set(block.values.keys()) == set(layers)

    @given(config=jvalue_block_config_st())
    @settings(max_examples=200)
    def test_sza_index_is_zero_for_single_block(self, config):
        """A single block always gets sza_index == 0."""
        species, layers, values = config

        stdout = build_jvalue_block(species, layers, values)
        result = parse_jvalues(stdout)

        assert len(result) == 1
        assert result[0].sza_index == 0
