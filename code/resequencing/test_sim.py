"""Verification tests for the Gate 1 simulator and oracle.

Run: python -I code/resequencing/test_sim.py
"""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim import draw_future, eval_seq, local_search, posterior_theta  # noqa: E402


def py_eval(seq, t0, T, Sv, opn, cls, w, wait):
    tot = 0.0
    for s in range(T.shape[0]):
        t, prev, late = t0, 0, 0.0
        for j in seq:
            t += T[s, prev, j]
            if wait and t < opn[j]:
                t = opn[j]
            late += w[j] * max(0.0, t - cls[j])
            t += Sv[s, j]
            prev = j
        tot += late
    return tot / T.shape[0]


def test_hand_example():
    # current node 0 -> stops 1, 2; deterministic; one scenario
    T = np.array([[[0, 10, 20], [10, 0, 5], [20, 5, 0]]], dtype=float)
    Sv = np.array([[0, 2, 3]], dtype=float)
    opn = np.array([-1e18, -1e18, -1e18])
    cls = np.array([1e18, 30.0, 12.0])
    w = np.array([0.0, 1.0, 1.0])
    # order 1,2: arrive 1 at 10 (ok), leave 12, arrive 2 at 17 -> late 5
    assert abs(eval_seq(np.array([1, 2]), 0.0, T, Sv, opn, cls, w, False) - 5.0) < 1e-9
    # order 2,1: arrive 2 at 20 -> late 8; leave 23, arrive 1 at 28 -> ok -> 8
    assert abs(eval_seq(np.array([2, 1]), 0.0, T, Sv, opn, cls, w, False) - 8.0) < 1e-9
    # waiting: open stop 1 at 15 -> arrive 10 wait to 15, leave 17, arrive 2 at 22 -> late 10
    opn2 = np.array([-1e18, 15.0, -1e18])
    assert abs(eval_seq(np.array([1, 2]), 0.0, T, Sv, opn2, cls, w, True) - 10.0) < 1e-9
    print("hand example ok")


def random_instance(rng, U, S, wait):
    N = U + 1
    xy = rng.uniform(0, 100, size=(N, 2))
    tb = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1)) * 3
    sb = rng.uniform(10, 60, size=N)
    T, Sv = draw_future(rng, tb, sb, 0.0, 0.02, 0.15, 0.25, 0.4, S)
    cls = rng.uniform(100, 900, size=N)
    opn = cls - rng.uniform(50, 400, size=N) if wait else np.full(N, -1e18)
    w = rng.integers(0, 2, size=N).astype(float)
    w[0] = 0.0
    return T, Sv, opn, cls, w


def test_eval_matches_python():
    rng = np.random.default_rng(1)
    for wait in (False, True):
        for _ in range(20):
            U = int(rng.integers(2, 12))
            T, Sv, opn, cls, w = random_instance(rng, U, 16, wait)
            seq = rng.permutation(np.arange(1, U + 1)).astype(np.int64)
            a = eval_seq(seq, 5.0, T, Sv, opn, cls, w, wait)
            b = py_eval(seq, 5.0, T, Sv, opn, cls, w, wait)
            assert abs(a - b) < 1e-6 * max(1, abs(b)), (a, b)
    print("eval_seq matches python reference")


def test_local_search_vs_bruteforce():
    rng = np.random.default_rng(7)
    gaps, hits, total = [], 0, 0
    for wait in (False, True):
        for _ in range(40):
            U = int(rng.integers(3, 8))
            T, Sv, opn, cls, w = random_instance(rng, U, 16, wait)
            cur = np.arange(1, U + 1, dtype=np.int64)
            start = eval_seq(cur, 0.0, T, Sv, opn, cls, w, wait)
            best_bf = min(eval_seq(np.array(p, dtype=np.int64), 0.0, T, Sv, opn, cls, w, wait)
                          for p in itertools.permutations(range(1, U + 1)))
            seq, val, _ = local_search(cur, 0.0, T, Sv, opn, cls, w, wait, 60)
            # reported value equals a re-evaluation of the returned sequence
            assert abs(val - eval_seq(seq, 0.0, T, Sv, opn, cls, w, wait)) < 1e-6 * max(1, val)
            assert sorted(seq.tolist()) == list(range(1, U + 1))
            assert val <= start + 1e-9
            total += 1
            if start > 0:
                gaps.append((val - best_bf) / max(start - best_bf, 1e-9) if start - best_bf > 1e-9 else 0.0)
            if val <= best_bf + 1e-6 * max(1, best_bf):
                hits += 1
    print(f"local search reached brute-force optimum in {hits}/{total} instances; "
          f"mean share of achievable improvement missed: {np.mean(gaps):.3f}")
    assert hits / total >= 0.8


def test_posterior():
    mu, v = posterior_theta(np.array([]), 0.15, 0.25)
    assert mu == 0.0 and abs(v - 0.0225) < 1e-12
    mu, v = posterior_theta(np.full(100, 0.2), 0.15, 0.25)
    assert 0.17 < mu < 0.2 and v < 0.0225
    print("posterior ok")


if __name__ == "__main__":
    test_hand_example()
    test_eval_matches_python()
    test_posterior()
    test_local_search_vs_bruteforce()
    print("ALL TESTS PASSED")
