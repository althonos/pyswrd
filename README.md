# 🐍🗡️ PySWRD [![Stars](https://img.shields.io/github/stars/althonos/pyswrd.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyswrd/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [SWORD](https://github.com/rvaser/sword) (Smith Waterman On Reduced Database), a method for fast database search.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyswrd/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyswrd/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyswrd?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pyswrd/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![PyPI](https://img.shields.io/pypi/v/pyswrd.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyswrd)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyswrd?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyswrd)
[![AUR](https://img.shields.io/aur/version/python-pyswrd?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyswrd)
[![Wheel](https://img.shields.io/pypi/wheel/pyswrd.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyswrd/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyswrd.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyswrd/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyswrd.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyswrd/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyswrd/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/pyswrd/)
[![Issues](https://img.shields.io/github/issues/althonos/pyswrd.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyswrd/issues)
[![Docs](https://img.shields.io/readthedocs/pyswrd/latest?style=flat-square&maxAge=600)](https://pyswrd.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyswrd/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyswrd?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyswrd)


## 🗺️ Overview

Searching a sequence inside a database of target sequences involves aligning
the sequence to all the targets to find the highest scoring ones, which has
a high computational cost. Several methods have been proposed over the years
that use a pre-filter to select. In BLAST[\[1\]](#ref1), k-mers are extracted 
from the query, and only targets containing high-scoring k-mers, with respect to 
the scoring matrix, are actually aligned.

[SWORD](https://github.com/rvaser/sword)[\[2\]](#ref2) proposes a pre-filter built 
on perfect hashing of short mismatching k-mers. The k-mers generated from the 
query sequence also include k-mers with mismatches to improve sensitivity. 
When a k-mer is found in a target sequence, SWORD computes the diagonal where it 
is located, similarly to FASTA[\[3\]](#ref3). Target sequences are then selected based on the
number of hits they have on the same diagonal. The pairwise alignment
is then handled by the platform-accelerated [Opal](https://github.com/Martinsos/opal)[\[4\]](#ref4)
library.

PySWRD is a [Python](https://python.org) module that provides bindings to
the heuristic filter part of [SWORD](https://github.com/rvaser/sword)
using [Cython](https://cython.org/). It implements a user-friendly, Pythonic
interface to build a heuristic filter, process a database in chunks, and
produce the indices of targets passing the filter. The resulting indices
can be used to filter a [PyOpal](https://github.com/althonos/pyopal) database,
using [Opal](https://github.com/Martinsos/opal) for pairwise alignment like
the original C++ implementation.

- **no binary dependency**: PySWRD is distributed as a Python package, so
  you can add it as a dependency to your project, and stop worrying about the
  SWORD binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you control, so you don't have to invoke the SWORD CLI using a sub-process
  and temporary files.
- **better portability**: Using only the heuristic filter of SWORD allows
  the code to be independent of the local CPU features, unlike SWORD and
  Opal which require SIMD. PySWRD delegates the SIMD compilation and
  dynamic dispatch to PyOpal to make the package easier to install. It
  also benefits from the wider platform support of PyOpal compared to
  the original Opal, featuring support for Windows and for Aarch64 CPUs.

## 🔧 Installing

PySWRD is available for all modern Python versions (3.6+).

It can be installed directly from [PyPI](https://pypi.org/project/pyopal/),
which hosts some pre-built x86-64 wheels for Linux, MacOS, and Windows, 
as well as the code required to compile from source with Cython:
```console
$ pip install pyswrd
```

<!-- Otherwise, PySWRD is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda pyswrd
``` -->

<!-- Check the [*install* page](https://pyswrd.readthedocs.io/en/stable/install.html)
of the documentation for other ways to install PySWRD on your machine. -->

## 💡 Example

PySWRD does not provide I/O, so the sequences to be used have to be loaded through
another library, such as [Biopython](https://biopython.org). PySWRD only requires 
the sequences to be available as Python strings:

```python
targets = [
    'MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPK', 
    'MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQTCASG', 
    'MASNTVSAQGGSNRPVRDFSNIQDVAQFLLFDPIWNEQPGSIVPWKMNRE', 
    'MYQAINPCPQSWYGSPQLEREIVCKMSGAPHYPNYYPVHPNALGGAWFDT', 
    'MARPLLGKTSSVRRRLESLSACSIFFFLRKFCQKMASLVFLNSPVYQMSN'
]
queries = [
    'MASNTVSAQGGSNRPVRDFSNIQDVAQFLLFDPIWNEQPG', 
    'MSFKVYDPIAELIATQFPTSNPDLQIINNDVLVVSPHKIT', 
    'MEQVPIKEMRLSDLRPNNKSIDTDLGGTKLVVIGKPGSGK'
]
```

Use the high-level `search` function, which wraps the internal classes in a single 
function to quickly run many-to-many searches in the event all your sequences are in 
memory. It expects the sequences as iterable of Python strings, and yields hits 
passing E-value and alignment thresholds:

```python
import pyswrd
for hit in pyswrd.search(queries, targets):
    print(hit.query_index, hit.target_index, hit.score, hit.evalue)
```

Different parameters can be passed to `pyswrd.search` and are passed to the 
SWORD filter and Opal alignment. For instance, to run SWORD in *fast* mode
instead of the default *sensitive* mode, and using the PAM70 matrix instead
of BLOSUM62, use:
```python
for hit in pyswrd.search(queries, targets, scorer_name="PAM70", score_threshold=0, kmer_length=5):
    print(hit.query_index, hit.target_index, hit.score, hit.evalue)
```

By default multithreading is supported, using one thread per CPU on the local
machine as reported by `os.cpu_count`, but it can be changed with the `threads` 
argument:
```python
for hit in pyswrd.search(queries, targets, threads=1):
    print(hit.query_index, hit.target_index, hit.score, hit.evalue)
```

You can also use the `pyswrd.HeuristicFilter` class directly if you wish to 
manage the data yourself, or if you want to use a different aligner.


<!-- ## 🧶 Thread-safety -->

<!-- ## ⏱️ Benchmarks -->


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue tracker](https://github.com/althonos/pyswrd/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.


### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyswrd/blob/main/CONTRIBUTING.md)
for more details.


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pyswrd/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).
SWORD was written by [Robert Vaser](https://github.com/rvaser) and is distributed under the terms of the
GPLv3 as well. See `vendor/sword/LICENSE` for more information. SWORD redistributes additional
libraries under the terms of the [MIT License](https://choosealicense.com/licenses/mit/).

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [SWORD authors](https://github.com/rvaser). It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*


## 📚 References

- <a id="ref1">\[1\]</a> Stephen F. Altschul, Warren Gish, Webb Miller, Eugene W. Myers, David J. Lipman. Basic local alignment search tool. J Mol Biol. 1990 Oct 5;215(3):403-10. [doi:10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2). [PMID:2231712](https://pubmed.ncbi.nlm.nih.gov/2231712).
- <a id="ref2">\[2\]</a> Robert Vaser, Dario Pavlović, Mile Šikić. SWORD—a highly efficient protein database search. Bioinformatics, Volume 32, Issue 17, September 2016, Pages i680–i684, [doi:10.1093/bioinformatics/btw445](https://doi.org/10.1093/bioinformatics/btw445).
- <a id="ref3">\[3\]</a> David J. Lipman, William R. Pearson. Rapid and sensitive protein similarity searches. Science. 1985 Mar 22;227(4693):1435-41. [doi:10.1126/science.2983426](https://doi.org/10.1126/science.2983426). [PMID:2983426](https://pubmed.ncbi.nlm.nih.gov/2983426).
- <a id="ref4">\[4\]</a> Korpar Matija, Martin Šošić, Dino Blažeka, Mile Šikić. SW#db: ‘GPU-Accelerated Exact Sequence Similarity Database Search’. PLoS One. 2015 Dec 31;10(12):e0145857. [doi:10.1371/journal.pone.0145857](https://doi.org/10.1371/journal.pone.0145857). [PMID:26719890](https://pubmed.ncbi.nlm.nih.gov/26719890). [PMC4699916](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4699916/).