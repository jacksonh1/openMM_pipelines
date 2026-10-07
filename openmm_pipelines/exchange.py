"""Replica-exchange mixing primitive, shared by the multistate backends and their analysis.

One place for the one computation every multistate caller needs — per-adjacent-rung exchange
acceptance — so the `read_mixing_statistics` PER-ITERATION-sum gotcha (see GOTCHAS.md) lives
once instead of in each backend's run report and the analysis module. Needs openmmtools (the
`[remd]` extra) but NOT mdtraj/matplotlib, so the produce side imports it without pulling the
analysis stack.
"""

from __future__ import annotations

import numpy as np


def neighbor_acceptance(reporter) -> list[float]:
    """Cumulative exchange acceptance for each adjacent rung pair i <-> i+1.

    `reporter` is an open `MultiStateReporter`. `read_mixing_statistics` returns PER-ITERATION
    matrices `(n_iter, n_states, n_states)` (NOT the documented `(n_states, n_states)`), so sum
    over iterations before indexing. Returns one value per neighbour pair (length n_states-1);
    a pair with no proposed swaps is reported as 0.0 (a broken ladder link).
    """
    n_accepted, n_proposed = reporter.read_mixing_statistics()
    n_accepted = np.asarray(n_accepted, dtype=float).sum(axis=0)
    n_proposed = np.asarray(n_proposed, dtype=float).sum(axis=0)
    n_states = n_accepted.shape[0]
    assert n_states >= 2, f"need >= 2 states for neighbour pairs, got {n_states}"
    return [
        float(n_accepted[i, i + 1] / n_proposed[i, i + 1]) if n_proposed[i, i + 1] > 0 else 0.0
        for i in range(n_states - 1)
    ]
