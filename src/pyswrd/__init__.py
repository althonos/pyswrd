# noqa: D104
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
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

from . import lib
from .lib import (
    __version__,
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

__doc__ = lib.__doc__

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

