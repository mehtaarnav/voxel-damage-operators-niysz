"""Validate cmlib/metrics.py against networks with analytically known answers.

  1. Series chain of k identical conductances c between the terminals:
        g_eff = c / k      (resistors in series)
        min-cut = c        (every single edge is a valid cut)
  2. m parallel chains:
        g_eff = m * c / k
        min-cut = m * c
  3. Deliberate bottleneck: two dense blobs joined by ONE weak edge.
        min-cut = weight of that edge
        g_eff  < min-cut   (series resistance inside the blobs)
  4. lambda_2 of an unweighted path graph P_n:
        lambda_2 = 2*(1 - cos(pi/n))     (exact, standard result)
"""
from __future__ import annotations

import os
import sys

import networkx as nx
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.metrics import (algebraic_connectivity, effective_conductance,  # noqa: E402
                           mincut_between_faces)

ok_all = True


def check(name, got, want, tol=1e-6):
    global ok_all
    good = abs(got - want) <= tol * max(1.0, abs(want))
    ok_all &= good
    print(f"  {name:46s} got={got:12.6f}  want={want:12.6f}  "
          f"{'PASS' if good else 'FAIL'}")


print("=" * 74)
print("TEST 1 — series chain")
for k in (2, 5, 10):
    c = 3.0
    G = nx.Graph()
    for i in range(k):
        G.add_edge(i, i + 1, cond=c)
    g = effective_conductance(G, [0], [k])
    mc, _ = mincut_between_faces(G, [0], [k])
    check(f"k={k:2d} g_eff", g, c / k)
    check(f"k={k:2d} min-cut", mc, c)

print("\n" + "=" * 74)
print("TEST 2 — m parallel chains of length k")
for m, k in ((2, 4), (3, 5)):
    c = 2.5
    G = nx.Graph()
    S, T = "S", "T"
    for j in range(m):
        prev = S
        for i in range(k - 1):
            node = (j, i)
            G.add_edge(prev, node, cond=c)
            prev = node
        G.add_edge(prev, T, cond=c)
    g = effective_conductance(G, [S], [T])
    mc, _ = mincut_between_faces(G, [S], [T])
    check(f"m={m},k={k} g_eff", g, m * c / k)
    check(f"m={m},k={k} min-cut", mc, m * c)

print("\n" + "=" * 74)
print("TEST 3 — single weak bottleneck between two dense blobs")
G = nx.Graph()
A = nx.complete_graph(range(0, 6))
B = nx.complete_graph(range(100, 106))
for u, v in A.edges:
    G.add_edge(u, v, cond=100.0)
for u, v in B.edges:
    G.add_edge(u, v, cond=100.0)
G.add_edge(5, 100, cond=0.7)          # the bottleneck
mc, note = mincut_between_faces(G, [0], [105])
g = effective_conductance(G, [0], [105])
check("min-cut equals the bottleneck edge", mc, 0.7)
print(f"  g_eff = {g:.6f}  (must be < min-cut = 0.7): "
      f"{'PASS' if g < 0.7 else 'FAIL'}")
ok_all &= (g < 0.7)

print("\n" + "=" * 74)
print("TEST 4 — lambda_2 of an unweighted path graph, exact formula")
for n in (4, 8, 16):
    G = nx.path_graph(n)
    nx.set_edge_attributes(G, 1.0, "cond")
    got = algebraic_connectivity(G, normalized=False)
    want = 2.0 * (1.0 - np.cos(np.pi / n))
    check(f"P_{n} lambda2", got, want, tol=1e-4)

print("\n" + "=" * 74)
print(f"METRICS VALIDATION: {'PASS' if ok_all else 'FAIL'}")
sys.exit(0 if ok_all else 1)
