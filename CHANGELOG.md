# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyswrd/compare/v0.2.0...HEAD


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
