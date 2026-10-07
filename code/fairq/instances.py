"""Instance generation: conflict graphs, rates, owners, colouring and independent sets."""
import numpy as np


def unit_disk_graph(N, density, seed):
    """N points uniform in a square of area N/density; edge when Euclidean distance <= 1."""
    rng = np.random.default_rng(seed)
    L = np.sqrt(N / density)
    pos = rng.uniform(0.0, L, size=(N, 2))
    edges = [(u, v) for u in range(N) for v in range(u + 1, N)
             if np.linalg.norm(pos[u] - pos[v]) <= 1.0]
    return edges, pos


def erdos_renyi_graph(N, p, seed):
    """G(N, p): each pair is an edge independently with probability p."""
    rng = np.random.default_rng(seed)
    return [(u, v) for u in range(N) for v in range(u + 1, N) if rng.random() < p], None


def random_regular_graph(N, d, seed):
    """Uniform random d-regular graph (networkx)."""
    import networkx as nx
    G = nx.random_regular_graph(d, N, seed=seed)
    return sorted(tuple(sorted(e)) for e in G.edges()), None


def neighbours(N, edges):
    nb = [set() for _ in range(N)]
    for u, v in edges:
        nb[u].add(v)
        nb[v].add(u)
    return [sorted(s) for s in nb]


def greedy_colouring(N, nb):
    """Largest-degree-first greedy proper colouring; returns the list of colour classes."""
    order = sorted(range(N), key=lambda v: -len(nb[v]))
    col = {}
    for v in order:
        used = {col[u] for u in nb[v] if u in col}
        c = 0
        while c in used:
            c += 1
        col[v] = c
    k = max(col.values()) + 1 if col else 0
    classes = [sorted(v for v in range(N) if col[v] == c) for c in range(k)]
    for cls in classes:                      # properness check
        for a in cls:
            assert not set(nb[a]) & set(cls)
    return classes


def independent_sets(N, nb):
    """All independent sets as sorted integer bitmasks (bit v set iff x_v = 1)."""
    nbmask = [sum(1 << u for u in nb[v]) for v in range(N)]
    out = []

    def rec(v, mask, forbidden):
        if v == N:
            out.append(mask)
            return
        rec(v + 1, mask, forbidden)
        if not (forbidden >> v) & 1:
            rec(v + 1, mask | (1 << v), forbidden | nbmask[v])

    rec(0, 0, 0)
    return np.array(sorted(out), dtype=np.int64)


FAMILIES = ("udg", "er", "reg3")


def make_instance(N, n_agents, seed, kind="udg", density=2.5, p=0.3):
    """Conflict graph of family `kind`, with rates r_v ~ U[0.5, 1.5] and balanced random owners."""
    if kind == "udg":
        edges, pos = unit_disk_graph(N, density, seed)
    elif kind == "er":
        edges, pos = erdos_renyi_graph(N, p, seed)
    elif kind == "reg3":
        edges, pos = random_regular_graph(N, 3, seed)
    else:
        raise ValueError(f"unknown family {kind!r}; expected one of {FAMILIES}")
    rng = np.random.default_rng(seed + 10_000)
    nb = neighbours(N, edges)
    owners = np.array([v % n_agents for v in rng.permutation(N)])
    rates = rng.uniform(0.5, 1.5, size=N)
    return dict(kind=kind, N=N, n=n_agents, edges=edges, pos=pos, nb=nb, owners=owners,
                rates=rates, colours=greedy_colouring(N, nb), F=independent_sets(N, nb))
