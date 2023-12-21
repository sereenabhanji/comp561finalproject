""" Graph class with Union-Find DS
and Kruskal's alg (finding max spanning tree) """

import random


def union(parent, rank, v1, v2):
    if rank[v1] < rank[v2]:
        parent[v1] = v2
    elif rank[v1] > rank[v2]:
        parent[v2] = v1
    else:
        parent[v2] = v1
        rank[v1] += 1


def print_graph(g):
    for e in g:
        print(str(e))


class Graph:
    edges = []

    def __init__(self, vertices):
        self.size = vertices

    def add_edge(self, e):
        self.edges.append(e)
        self.edges.sort(key=lambda tup: tup[2])  # sorted ascending, will pop to get largest

    def kruskal(self):
        mst = []
        s = 0

        parent = []
        rank = []
        for v in range(self.size):
            parent.append(v)
            rank.append(0)

        while len(self.edges) > 0:
            curr = self.edges.pop()
            s += 1
            u1 = self.find(parent, curr[0])
            u2 = self.find(parent, curr[1])

            if u1 != u2:
                mst.append(curr)
                union(parent, rank, u1, u2)

        print_graph(mst)
        return mst

    def find(self, parent, i):
        if parent[i] != i:
            parent[i] = self.find(parent, parent[i])
        return parent[i]


if __name__ == '__main__':
    g = Graph(5)
    for i in range(4):
        w = random.randint(-5, 12)
        edge = (i, i + 1, w)
        g.add_edge(edge)
        g.kruskal()
