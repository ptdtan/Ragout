import networkx as nx
from collections import namedtuple
import os

from permutation import *
from debug import DebugConfig, write_dot
import phylogeny as phylo


Connection = namedtuple("Connection", ["start", "end"])


class BreakpointGraph:
    def __init__(self):
        self.bp_graph = nx.MultiGraph()
        self.targets = []
        self.references = []
        self.known_adjacencies = {}


    def build_from(self, perm_container, circular):
        for perm in perm_container.ref_perms_filtered:
            if perm.ref_id not in self.references:
                self.references.append(perm.ref_id)

            prev = None
            for block in perm.iter_blocks(circular):
                if not prev:
                    prev = block
                    continue

                left_block = prev
                right_block = block
                self.bp_graph.add_node(-left_block)
                self.bp_graph.add_node(right_block)
                self.bp_graph.add_edge(-left_block, right_block, ref_id=perm.ref_id)

                prev = block

        for perm in perm_container.target_perms_filtered:
            if perm.ref_id not in self.targets:
                self.targets.append(perm.ref_id)

            prev = None
            for block in perm.iter_blocks(False):
                if not prev:
                    prev = block
                    continue

                self.known_adjacencies[-prev] = block
                self.known_adjacencies[block] = -prev
                prev = block


    def find_adjacencies(self, phylogeny):
        adjacencies = {}
        subgraphs = nx.connected_component_subgraphs(self.bp_graph)
        for comp_id, subgraph in enumerate(subgraphs):
            if len(subgraph) < 2:
                continue

            if len(subgraph) == 2:
                node_1, node_2 = subgraph.nodes()
                adjacencies[-node_1] = Connection(-node_1, node_2)
                adjacencies[-node_2] = Connection(-node_2, node_1)
                continue

            weighted_graph = self.make_weighted(subgraph, phylogeny)
            chosen_edges = split_graph(weighted_graph)

            for edge in chosen_edges:
                adjacencies[-edge[0]] = Connection(-edge[0], edge[1])
                adjacencies[-edge[1]] = Connection(-edge[1], edge[0])

            if DebugConfig.get_writer().debugging:
                debug_dir = DebugConfig.get_writer().debug_dir
                debug_draw_component(comp_id, weighted_graph, subgraph, debug_dir)

        return adjacencies


    def make_weighted(self, graph, phylogeny):
        assert len(graph) > 2
        g = nx.Graph()
        g.add_nodes_from(graph.nodes())
        target_id = self.targets[0]

        for node in graph.nodes():
            adjacencies = {}
            for neighbor in graph.neighbors(node):
                for edge in graph[node][neighbor].values():
                    adjacencies[edge["ref_id"]] = neighbor

            for ref_id in self.references:
                if ref_id not in adjacencies:
                    adjacencies[ref_id] = None  #"void" state in paper

            for neighbor in graph.neighbors(node):
                break_weight = 0.0
                if not (self.known_adjacencies.get(node, None) == neighbor):
                    adjacencies[target_id] = neighbor
                    break_weight = phylogeny.estimate_tree(adjacencies)

                update_edge(g, node, neighbor, break_weight)

        return g


def split_graph(graph):
    for v1, v2 in graph.edges_iter():
        graph[v1][v2]["weight"] = -graph[v1][v2]["weight"] #want minimum weight

    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = []
    for v1, v2 in edges.iteritems():
        if not (v2, v1) in unique_edges:
            unique_edges.append((v1, v2))

    return unique_edges


def update_edge(graph, v1, v2, weight):
    if not graph.has_edge(v1, v2):
        graph.add_edge(v1, v2, weight=weight)
    else:
        graph[v1][v2]["weight"] += weight


def debug_draw_component(comp_id, weighted_graph, breakpoint_graph, debug_dir):
    if len(breakpoint_graph) == 2:
        return

    for e in weighted_graph.edges_iter():
        weighted_graph[e[0]][e[1]]["label"] = ("{0:7.4f}"
                                    .format(weighted_graph[e[0]][e[1]]["weight"]))
    bg_out = os.path.join(debug_dir, "comp{0}-bg.dot".format(comp_id))
    weighted_out = os.path.join(debug_dir, "comp{0}-weighted.dot".format(comp_id))
    write_dot(breakpoint_graph, open(bg_out, "w"))
    write_dot(weighted_graph, open(weighted_out, "w"))
