from collections import defaultdict


class IntervalNode:
    def __init__(self):
        self.values = set()
        self.neighbors = set()

    def __str__(self):
        value_string = "values: " + str(self.values)
        neighbor_string = "neighbors: " + str(self.neighbors)

        return value_string + '\n' + neighbor_string


class IntervalGraph:
    def __init__(self, interval_dict: dict):
        self.graph = defaultdict(lambda: IntervalNode)
        self.construct_from_intervals(interval_dict)

    def construct_from_intervals(self, interval_dict: dict):
        tuples = [(x[0],x[1],y) for x,y in interval_dict.items()]

        # Sort the intervals based on their start points
        sorted_intervals = sorted(tuples, key=lambda x: x[0])

        active_intervals = []

        # Sweep through the sorted intervals
        for start, stop, value in sorted_intervals:
            interval = (start,stop)

            if interval not in self.graph:
                self.graph[interval] = IntervalNode()
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


def test():
    intervals = {
        (1,3):"a",
        (2,3):"b",
        (2,4):"c",
        (4,5):"c",
        (6,7):"c"
    }

    graph = IntervalGraph(intervals)

    for interval,node in graph.graph.items():
        print(interval)
        print(node)
        print()


if __name__ == "__main__":
    test()

