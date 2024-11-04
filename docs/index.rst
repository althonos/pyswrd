PySWRD |Stars|
==============

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyswrd.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyswrd/stargazers

`Cython <https://cython.org/>`_ *bindings and Python interface to* `SWORD <https://github.com/rvaser/sword>`_ *(Smith Waterman On Reduced Database), a method for fast database search.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyswrd/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyswrd/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyswrd?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyswrd/

.. |PyPI| image:: https://img.shields.io/pypi/v/pyswrd.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyswrd

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyswrd?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyswrd

.. |AUR| image:: https://img.shields.io/aur/version/python-pyswrd?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyswrd

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyswrd?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyswrd/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyswrd.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyswrd/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyswrd.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyswrd/#files

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyswrd/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/pyswrd/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyswrd.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyswrd/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pyswrd?style=flat-square&maxAge=3600
   :target: http://pyswrd.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyswrd/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyswrd?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pyswrd


Overview
--------

Searching a sequence inside a database of target sequences involves aligning
the sequence to all the targets to find the highest scoring ones, which has
a high computational cost. Several methods have been proposed over the years
that use a pre-filter to select. In BLAST, k-mers are extracted 
from the query, and only targets containing high-scoring k-mers, with respect to 
the scoring matrix, are actually aligned.

`SWORD <https://github.com/rvaser/sword>`_ proposes a pre-filter built 
on perfect hashing of short mismatching k-mers. The k-mers generated from the 
query sequence also include k-mers with mismatches to improve sensitivity. 
When a k-mer is found in a target sequence, SWORD computes the diagonal where it 
is located, similarly to FASTA. Target sequences are then selected based on the
number of hits they have on the same diagonal. The pairwise alignment
is then handled by the platform-accelerated 
`Opal <https://github.com/Martinsos/opal>`_ library.

PySWRD is a `Python <https://python.org>`_ module that provides bindings to
the heuristic filter part of `SWORD <https://github.com/rvaser/sword>`_
using `Cython <https://cython.org/>`_. It implements a user-friendly, Pythonic
interface to build a heuristic filter, process a database in chunks, and
produce the indices of targets passing the filter. The resulting indices
can be used to filter a `PyOpal <https://github.com/althonos/pyopal>`_ database,
using `Opal <https://github.com/Martinsos/opal>`_ for pairwise alignment like
the original C++ implementation.


Setup
-----

Run ``pip install pyswrd`` in a shell to download the latest release from PyPi, 
or have a look at the :doc:`Installation page <guide/install>` to find other ways 
to install ``pyswrd``.


Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst
   

License
-------

This library is provided under the `GNU General Public License 3.0 or later <https://choosealicense.com/licenses/gpl-3.0/>`_.
SWORD was developed by `Robert Vaser <https://github.com/rvaser>`_ and 
`Dario Pavlovic <https://github.com/darxsys>`_ under the terms of the 
`GNU General Public License 3.0 or later <https://choosealicense.com/licenses/gpl-3.0/>`_
as well. See the :doc:`Copyright Notice <guide/copyright>` section
for the full license.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `SWORD`_ *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
