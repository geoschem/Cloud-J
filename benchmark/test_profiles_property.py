"""Property-based tests for benchmark.profiles module.

Validates: Requirements 2.3, 2.4
"""

from hypothesis import given, settings
from hypothesis import strategies as st

from benchmark.profiles import build_scenario_matrix
from benchmark.types import ProfileCondition, Scenario

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

# Strategy for generating ProfileCondition objects with unique names
profile_condition_st = st.builds(
    ProfileCondition,
    name=st.text(min_size=1, max_size=20, alphabet=st.characters(whitelist_categories=("L", "N", "Pd"))),
    description=st.text(min_size=0, max_size=50),
    cloud_scale=st.floats(min_value=0.0, max_value=10.0, allow_nan=False, allow_infinity=False),
    aerosol_scale=st.floats(min_value=0.0, max_value=10.0, allow_nan=False, allow_infinity=False),
)

# Strategy for non-empty lists of profiles with unique names
profiles_st = st.lists(
    profile_condition_st,
    min_size=1,
    max_size=10,
    unique_by=lambda p: p.name,
)

# Strategy for SZA values (valid solar zenith angles)
sza_st = st.floats(min_value=0.0, max_value=90.0, allow_nan=False, allow_infinity=False)

# Strategy for non-empty lists of unique SZA values
szas_st = st.lists(
    sza_st,
    min_size=1,
    max_size=10,
    unique=True,
)


# ---------------------------------------------------------------------------
# Property 1: Scenario matrix is the Cartesian product
# ---------------------------------------------------------------------------


class TestScenarioMatrixProperty:
    """**Validates: Requirements 2.3, 2.4**"""

    @given(profiles=profiles_st, szas=szas_st)
    @settings(max_examples=200)
    def test_length_equals_product(self, profiles, szas):
        """The result length equals len(profiles) * len(szas)."""
        result = build_scenario_matrix(profiles, szas)
        assert len(result) == len(profiles) * len(szas)

    @given(profiles=profiles_st, szas=szas_st)
    @settings(max_examples=200)
    def test_every_combination_appears_exactly_once(self, profiles, szas):
        """Every (profile, sza) combination appears exactly once."""
        result = build_scenario_matrix(profiles, szas)

        # Build expected set of (profile_name, sza) tuples
        expected = {(p.name, s) for p in profiles for s in szas}

        # Build actual set of (profile_name, sza) tuples
        actual = [(scenario.profile.name, scenario.sza) for scenario in result]

        # Check every expected combination is present
        actual_set = set(actual)
        assert actual_set == expected

    @given(profiles=profiles_st, szas=szas_st)
    @settings(max_examples=200)
    def test_no_duplicates(self, profiles, szas):
        """No duplicate scenarios exist in the result."""
        result = build_scenario_matrix(profiles, szas)

        # Extract (profile_name, sza) tuples
        pairs = [(scenario.profile.name, scenario.sza) for scenario in result]

        # No duplicates: length of list equals length of set
        assert len(pairs) == len(set(pairs))

    @given(profiles=profiles_st, szas=szas_st)
    @settings(max_examples=200)
    def test_all_results_are_scenario_objects(self, profiles, szas):
        """Every element in the result is a Scenario with correct types."""
        result = build_scenario_matrix(profiles, szas)
        for scenario in result:
            assert isinstance(scenario, Scenario)
            assert isinstance(scenario.profile, ProfileCondition)
            assert isinstance(scenario.sza, float)
