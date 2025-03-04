# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyswrd/compare/v0.3.0-post1...HEAD


## [v0.3.0-post1] - 2025-03-04
[v0.3.0-post1]: https://github.com/althonos/pyswrd/compare/v0.3.0...v0.3.0-post1

### Fixed
- Extra key in `pyproject.toml` causing build issues with version `0.11.0` of `scikit-build-core`.


## [v0.3.0] - 2024-11-04
[v0.3.0]: https://github.com/althonos/pyswrd/compare/v0.2.0...v0.3.0

### Added
- Support for Python 3.13.
- Type hints for `pyswrd` using external stub file.

### Changed
- Bump required `scoring-matrices` dependency to `v0.3.0`.
- Bump required `pyopal` depdency to `v0.7.0`.
- Reorganize project to build with CMake and `scikit-build-core`.
- Update documentation to use the PyData theme.

### Fixed
- Signature of Cython classes constructors not displaying in documentation.
- Large test data being included in source distribution.
- Incorrect license tag in Python files.

### Removed
- Support for Python 3.6.


## [v0.2.0] - 2024-05-15
[v0.2.0]: https://github.com/althonos/pyswrd/compare/v0.1.1...v0.2.0

### Changed
- Update `pyopal` dependency to `v0.6.0`.
- Use `scoring-matrices` package to expose the scoring matrix of `pyswrd.Scorer` objects.
- Inline some SWORD functions to improve performance in some critical sections.


## [v0.1.1] - 2024-02-28
[v0.1.1]: https://github.com/althonos/pyswrd/compare/v0.1.0...v0.1.1

### Added
- Support for passing a running `ThreadPool` to `pyswrd.search`.

### Changed
- Update `pyopal` dependency to `v0.5.2`.


## [v0.1.0] - 2024-01-20
[v0.1.0]: https://github.com/althonos/pyswrd/compare/1832ee2...v0.1.0

Initial release.
