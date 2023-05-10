from collections import defaultdict
import re

from matplotlib import pyplot
from networkx import DiGraph
import networkx


def plot_graph(graph, color_map, output_path):
    for layer, nodes in enumerate(networkx.topological_generations(graph)):
        # `multipartite_layout` expects the layer as a node attribute, so add the
        # numeric layer value as a node attribute
        for node in nodes:
            graph.nodes[node]["layer"] = layer

    pos = networkx.multipartite_layout(graph, subset_key="layer")

    n_ref = graph.number_of_nodes()
    width = 1 + (n_ref*0.35)
    height = 1 + (n_ref*0.07)
    f = pyplot.figure(figsize=(width, height), dpi=200)
    a = pyplot.axes()

    networkx.draw(
        graph,
        pos,
        connectionstyle="arc3,rad=-0.35",
        style=':',
        node_color=color_map,
        node_size=100,
        font_size=8,
        width=0.6,
        arrowsize=3,
        with_labels=True)

    for edge in graph.edges(data='weight'):
        print(edge)

        networkx.draw_networkx_edges(
            graph,
            pos,
            edgelist=[edge],
            width=edge[2],
            connectionstyle="arc3,rad=-0.35",
            arrowstyle='-',
            alpha=0.6,
            edge_color="#000000"
        )

    ylim = list(a.get_ylim())
    ylim[0] -= 0.5
    a.set_ylim(ylim)

    pyplot.savefig(output_path, dpi=200)


def main():
    # gfa_path = "/home/ryan/data/test_hapslap/output/chr20:47474020-47477018_no_empty.gfa"
    # gaf_path = "/home/ryan/data/test_hapslap/output/test.gaf"
    # csv_path = "/home/ryan/data/test_hapslap/output/chr20:47474020-47477018.csv"

    gfa_path = "/home/ryan/data/test_hapslap/output/chr20:54974920-54977307_no_empty.gfa"
    gaf_path = "/home/ryan/data/test_hapslap/output/test.gaf"
    csv_path = "/home/ryan/data/test_hapslap/output/chr20:54974920-54977307.csv"

    graph = DiGraph()
    max_id = -1

    with open(gfa_path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith("S"):
                id = line.strip().split('\t')[1]
                graph.add_node(id)

                if int(id) > max_id:
                    max_id = int(id)

        file.seek(0)

        for l,line in enumerate(file):
            if line.startswith("L"):
                tokens = line.strip().split('\t')
                a = tokens[1]
                b = tokens[3]
                a_reversal = tokens[2] == '-'
                b_reversal = tokens[4] == '-'

                print(a,b)

                if a_reversal == False and b_reversal == False:
                    graph.add_edge(a,b,weight=0)
                else:
                    exit("ERROR: non forward edge found in GFA line %d: %s" % (l,line))

    color_map = dict()
    color_index = None

    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split(',')

            print(tokens)
            if l == 0:
                color_index = tokens.index("color")
                continue

            id = tokens[0]
            color = tokens[color_index]
            color_map[id] = color

    color_list = list()

    for n in graph.nodes:
        if n in color_map:
            color_list.append(color_map[n])

    print(color_list)

    png_path = gfa_path.replace(".gfa", "_aligned.png")

    max_weight = 0.0
    with open(gaf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')

            path = re.split('[><]', tokens[5][1:])

            print(path)

            for i in range(len(path) - 1):
                a = path[i]
                b = path[i+1]

                print(a,b)

                if graph.has_edge(a,b):
                    print("found edge: a,b")
                    graph[a][b]["weight"] += 1
                    if graph[a][b]["weight"] > max_weight:
                        max_weight = graph[a][b]["weight"]

                if graph.has_edge(b,a):
                    print("found edge: b,a")
                    graph[b][a]["weight"] += 1
                    if graph[b][a]["weight"] > max_weight:
                        max_weight = graph[b][a]["weight"]

    max_width = 5.0

    for edge in graph.edges:
        graph[edge[0]][edge[1]]["weight"] /= max_weight/max_width

    plot_graph(graph=graph, color_map=color_list, output_path=png_path)


if __name__ == "__main__":
    main()

