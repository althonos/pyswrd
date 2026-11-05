import contextlib
import gzip
import io
import os

try:
    try:
        from importlib.resources import files
    except ImportError:
        from importlib_resources import files  # type: ignore
except ImportError:
    files = None  # type: ignore

from ..fasta import parse, zopen


def exists(name):
    return files(__name__).joinpath(name).is_file()


@contextlib.contextmanager
def load(name, mode="rb"):
    with files(__name__).joinpath(name).open(mode="rb") as src:
        with zopen(src, mode=mode) as reader:
            yield reader


def load_record(name):
    with load(name, mode="r") as f:
        return next(parse(f))


def load_records(name):
    with load(name, mode="r") as f:
        return list(parse(f))


def load_text(name):
    with load(name, mode="r") as f:
        return f.read()
