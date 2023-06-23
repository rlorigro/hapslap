from collections import defaultdict
from queue import SimpleQueue


class IntervalNode:
    def __init__(self):
        self.values = set()
        self.neighbors = set()

    def __str__(self):
        value_string = "values: " + str(self.values)
        neighbor_string = "neighbors: " + str(self.neighbors)

        return value_string + '\n' + neighbor_string


class IntervalGraph:
    def __init__(self, intervals: list):
        self.graph = dict()
        self.construct_from_intervals(intervals)

    def construct_from_intervals(self, intervals: list):
        # Sort the intervals based on their start points
        sorted_intervals = sorted(intervals, key=lambda x: x[0])

        active_intervals = []

        # Sweep through the sorted intervals
        for start, stop, value in sorted_intervals:
            if not (type(start) == int or type(start == float)):
                raise Exception("ERROR: start index for tuple is not integer or float: " + str(tuple))

            if not (type(stop) == int or type(stop == float)):
                raise Exception("ERROR: stop index for tuple is not integer or float: " + str(tuple))

            interval = (start,stop)

            if interval not in self.graph:
                self.graph[interval] = IntervalNode()

            # If the value is set or list, then add each item independently, but NOT if tuple (is hashable)
            if type(value) == set or type(value) == list:
                for v in value:
                    self.graph[interval].values.add(v)
            else:
                self.graph[interval].values.add(value)

            # Remove inactive intervals
            # TODO: this would be better performed with a heap DS, sorting by end, or some kind of tandem iterator
            active_intervals = [i for i in active_intervals if i[1] > start]

            # Add current interval to active intervals
            active_intervals.append(interval)

            # Connect current interval to all active intervals
            for other_interval in active_intervals:
                if other_interval != interval:
                    self.graph[other_interval].neighbors.add(interval)
                    self.graph[interval].neighbors.add(other_interval)

        return self.graph

    def get_connected_components(self):
        unvisited = set(self.graph.keys())
        components = list()

        while len(unvisited) > 0:
            components.append(set())

            q = SimpleQueue()
            q.put(unvisited.pop())

            while not q.empty():
                n = q.get()
                components[-1].add(n)

                for other in self.graph[n].neighbors:
                    if other in unvisited:
                        q.put(other)
                        unvisited.remove(other)

        return components

    def __str__(self):
        s = ""
        for interval,node in self.graph.items():
            s += str(interval)
            s += '\n'
            s += str(node)
            s += '\n'
            s += '\n'

        return s


def test():
    intervals = [
        (1,3,"a"),
        (1,3,"a2"),
        (2,3,"b"),
        (2,3,"b"),
        (2,4,"c"),
        (4,5,"d"),
        (6,7,"e"),
        (6,7,{"f","g"})
    ]

    graph = IntervalGraph(intervals)

    print(graph)

    components = graph.get_connected_components()

    for c,component in enumerate(components):
        print("---")
        print(c)
        print(component)


if __name__ == "__main__":
    test()
