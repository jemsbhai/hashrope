"""Shared fixtures for hashrope tests."""

import pytest
from hashrope import PolynomialHash, MERSENNE_61


@pytest.fixture
def h():
    """Default PolynomialHash instance (base=131, p=2^61-1)."""
    return PolynomialHash(prime=MERSENNE_61, base=131)


@pytest.fixture
def h257():
    """PolynomialHash with base=257 for multi-hash tests."""
    return PolynomialHash(prime=MERSENNE_61, base=257)
