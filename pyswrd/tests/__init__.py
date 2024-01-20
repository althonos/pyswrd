# noqa: D104

from . import (
    test_doctest,
    test_filter,
    test_search,
    test_sequences,
)


def load_tests(loader, suite, pattern):
    test_doctest.load_tests(loader, suite, pattern)
    suite.addTests(loader.loadTestsFromModule(test_filter))
    suite.addTests(loader.loadTestsFromModule(test_search))
    suite.addTests(loader.loadTestsFromModule(test_sequences))
    return suite
