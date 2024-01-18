# noqa: D104
from ._version import __version__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"

from pyopal import Aligner

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


def search(
    queries,
    targets,
    *,
    gap_open = 10,
    gap_extend = 1,
    scorer_name = "BLOSUM62",
    kmer_length = 3,
    max_candidates = 30000,
    score_threshold = 13,
    max_alignments = 10,
    max_evalue=10.0,
    algorithm = "sw",
    threads = 1,
):
    """Run a many-to-many search of query sequences to target sequences.

    This function is a high-level wrapper around the different classes of
    the `pyswrd` library to support fast searches when all sequences are
    in memory.

    Arguments:
        queries (`~pyswrd.Sequences`, or iterable of `str`): The sequences
            to query the target sequences with.
        targets (`~pyswrd.Sequences` or iterable of `str`): The sequences
            to be queries with the query sequences.
        gap_open (`int`): The penalty for opening a gap in each alignment.
        gap_extend (`int`): The penalty for extending a gap in each
            alignment.
        scorer_name (`str`): The name of the scoring matrix to use
            for scoring each alignment. See `~pyswrd.Scorer` for the
            list of supported names.
        kmer_length (`int`): The length of the k-mers to use in the
            SWORD heuristic filter.
        max_candidates (`int`): The maximum number of candidates to
            retain in the heuristic filter.
        max_evalue (`float`): The E-value threshold above which to
            discard sequences before alignment.
        algorithm (`str`): The algorithm to use to perform pairwise
            alignment. See `pyopal.Database.search` for more
            information.
        threads (`int`): The number of threads to use to run the
            pre-filter. When 0 is given, use all threads on the
            local machine as reported by `os.cpu_count`.

    Yields:
        `~pyswrd.Hit`: Hit objects for each hit passing the threshold
        parameters. Hits are grouped by query index, and sorted by
        E-value.

    Example:
        >>> queries = ["MAGFLKVVQLLAKYGSKAVQWAWANKGKILDWLNAGQAIDWVVSKIKQILGIK"]
        >>> targets = ([
        ...     "MESILDLQELETSEEESALMAASTVSNNC",
        ...     "MKKAVIVENKGCATCSIGAACLVDGPIPDFEIAGATGLFGLWG",
        ...     "MAGFLKVVQILAKYGSKAVQWAWANKGKILDWINAGQAIDWVVEKIKQILGIK",
        ...     "MTQIKVPTALIASVHGEGQHLFEPMAARCTCTTIISSSSTF",
        ... ])
        >>> for hit in pyswrd.search(queries, targets):
        ...     cigar = hit.result.cigar()
        ...     print(f"target={hit.target_index} score={hit.score} evalue={hit.evalue:.1g} cigar={cigar}")
        target=2 score=268 evalue=1e-33 cigar=53M

    """
    query_db  = queries if isinstance(queries, Sequences) else Sequences(queries)
    target_db = targets if isinstance(targets, Sequences) else Sequences(targets)
    target_lengths = target_db.lengths

    scorer  = Scorer(name=scorer_name, gap_open=gap_open, gap_extend=gap_extend)
    score_matrix = scorer.score_matrix
    hfilter = HeuristicFilter(query_db, kmer_length=kmer_length, max_candidates=max_candidates, score_threshold=score_threshold, scorer=scorer, threads=threads)
    aligner = Aligner(score_matrix, gap_open=gap_open, gap_extend=gap_extend)

    filter_result = hfilter.score(target_db).finish()
    evalue = EValue(filter_result.database_length, scorer)

    for query_index, query in enumerate(query_db):
        # get length of query
        query_length = len(query)
        # extract candidates and align them in scoring mode only
        target_indices = filter_result.indices[query_index]
        sub_db = target_db.extract(target_indices)
        score_results = aligner.align(query, sub_db, algorithm=algorithm, mode="score")
        # extract indices with E-value under threshold
        target_evalues = []
        for result, target_index in zip(score_results, target_indices):
            target_length = target_lengths[target_index]
            target_evalue = evalue.calculate(result.score, query_length, target_length)
            if target_evalue <= max_evalue:
                target_evalues.append((target_index, target_evalue))
        # get only `max_alignments` alignments per query, smallest e-values first
        target_evalues.sort(key=lambda x: x[1])
        target_indices = [x[0] for x in target_evalues]
        target_evalues = target_evalues[:max_alignments]
        # align selected sequences
        sub_db = target_db.extract(target_indices)
        ali_results = aligner.align(query, sub_db, algorithm=algorithm, mode="full")
        # return hits for aligned sequences
        for (target_index, target_evalue), target_result in zip(target_evalues, ali_results):
            yield Hit(
                query_index=query_index,
                target_index=target_index,
                evalue=target_evalue,
                result=target_result
            )