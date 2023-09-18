import argparse
import networkx
from copy import copy


def path_recursion(graph, id, path=None):
    if path is None:
        path = list()

    path.append(id)
    out_edges = graph.out_edges(id)

    if len(out_edges) == 0:
        yield tuple(path)
    else:
        for edge in out_edges:
            yield from path_recursion(graph=graph, id=edge[1], path=copy(path))


def enumerate_paths(graph):
    # Get start node
    start_id = next(networkx.topological_sort(graph))

    print("Starting path recursion from %s" % start_id)

    paths = [p for p in path_recursion(graph=graph, id=start_id)]

    return paths


def read_gfa_as_digraph(gfa_path, graph):
    with open(gfa_path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith("S"):
                data = line.strip().split('\t')
                id = data[1]
                sequence = data[2]

                print(id, sequence if len(sequence) <= 100 else sequence[:100])

                graph.add_node(id, sequence=sequence)

        file.seek(0)

        for l,line in enumerate(file):
            if line.startswith("L"):
                tokens = line.strip().split('\t')
                a = tokens[1]
                b = tokens[3]
                a_reversal = tokens[2] == '-'
                b_reversal = tokens[4] == '-'

                print(a,a_reversal,b,b_reversal)

                if a_reversal == False and b_reversal == False:
                    graph.add_edge(a,b)
                else:
                    exit("ERROR: non forward edge found in GFA line %d: %s" % (l,line))

    return graph


def path_to_sequence(graph, path):
    sequences = list()
    for id in path:
        sequences.append(graph.nodes[id]["sequence"])

    return ''.join(sequences)


def validate_graph_subset(path_a, path_b):
    graph_a = networkx.DiGraph()
    graph_b = networkx.DiGraph()

    graph_a = read_gfa_as_digraph(gfa_path=path_a, graph=graph_a)
    graph_b = read_gfa_as_digraph(gfa_path=path_b, graph=graph_b)

    paths_a = enumerate_paths(graph_a)
    paths_b = enumerate_paths(graph_b)

    path_seqs_a = {path_to_sequence(graph=graph_a, path=p):p for p in paths_a}
    path_seqs_b = {path_to_sequence(graph=graph_b, path=p):p for p in paths_b}

    seq_set_a = path_seqs_a.keys()
    seq_set_b = path_seqs_b.keys()

    failing_paths = set()
    for s in seq_set_a:
        path = path_seqs_a[s]
        print("validating path from graph A:", path)

        if s not in seq_set_b:
            failing_paths.add(path)

    if len(failing_paths) != 0:
        for p in failing_paths:
            print(p)

        raise Exception("FAIL: not all paths in graph A contained in graph B")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a",
        required=True,
        type=str,
        help="The graph (GFA) that is expected to be a subset of graph B"
    )

    parser.add_argument(
        "-b",
        required=True,
        type=str,
        help="The graph (GFA) that is expected to be a superset of graph A"
    )

    args = parser.parse_args()

    validate_graph_subset(path_a=args.a, path_b=args.b)
