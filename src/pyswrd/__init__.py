# noqa: D104
from ._version import __version__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__all__ = [
    "HeuristicFilter",
    "KmerGenerator",
    "Scorer",
    "EValue",
    "Sequences",
    "FilterScore",
    "FilterResult",
    "Hit",
    "search",
]

from . import _sword
from ._sword import (
    HeuristicFilter,
    KmerGenerator,
    Scorer,
    EValue,
    Sequences,
    FilterScore,
    FilterResult,
    Hit,
    search
)

__doc__ = _sword.__doc__

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pyswrd.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )

