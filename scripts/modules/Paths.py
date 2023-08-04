from copy import copy


class Paths:
    def __init__(self):
        self.path_to_id = dict()
        self.id_to_path = list()
        self.id_to_weight = list()
        self.n = 0

    def add_path(self, path: tuple[int], weight) -> int:
        p = None

        if path not in self.path_to_id:
            p = len(self.id_to_path)

            self.id_to_path.append(path)
            self.id_to_weight.append(weight)

            self.path_to_id[path] = p
            self.n += 1
        else:
            # Don't duplicate values if it already exists
            p = self.path_to_id[path]

        return p

    '''
    Increment an existing path weight (+= 'weight') or create a new path and initialize with 'weight' 
    '''
    def increment_weight(self, path, weight):
        if self.has_path(path):
            self.id_to_weight[self.get_path_id(path)] += weight
        else:
            self.add_path(path, weight)

    def get_path(self, path_id: int):
        if self.has_path(path_id):
            return self.id_to_path[path_id]
        else:
            raise Exception("ERROR: path ID %d not in paths" % path_id)

    def has_path(self, path):
        has_path = None

        if type(path) == tuple:
            has_path = (path in self.path_to_id)
        elif type(path) == str:
            has_path = (self.get_path_tuple(path) in self.path_to_id)
        elif type(path) == int:
            has_path = (path < self.n)
        else:
            raise Exception("ERROR: path of type %s not str (a_b_c), tuple[int] (a,b,c), or int (ID): %s" % (str(type(path)), str(path)))

        return has_path

    def get_path_id(self, path) -> int:
        p = None

        if type(path) == tuple:
            p = self.path_to_id[path]
        elif type(path) == str:
            p = self.path_to_id[self.get_path_tuple(path)]
        elif type(path) == int:
            p = path
        else:
            raise Exception("ERROR: path of type %s not str (a_b_c), tuple[int] (a,b,c), or int (ID): %s" % (str(type(path)), str(path)))

        return p

    def get_path_name(self, path) -> str:
        p = None

        if type(path) == tuple:
            p = '_'.join(list(map(str,path))) if len(path) > 1 else str(path[0])
        elif type(path) == str:
            p = path
        elif type(path) == int:
            p = '_'.join(list(map(str,self.id_to_path[path])))
        else:
            raise Exception("ERROR: path of type %s not str (a_b_c), tuple[int] (a,b,c), or int (ID): %s" % (str(type(path)), str(path)))

        return p

    def get_path_weight(self, path):
        return self.id_to_weight[self.get_path_id(path)]

    @staticmethod
    def get_path_tuple(path_name: str) -> tuple[int]:
        return tuple(map(int,path_name.split('_')))

    def ids(self):
        for p in range(self.n):
            yield p

    def __iter__(self):
        for p in range(self.n):
            yield p, self.id_to_path[p], self.id_to_weight[p]

    def __len__(self):
        return self.n


def path_recursion(graph, alleles, id, path_sequence=None):
    if path_sequence is None:
        path_sequence = list()

    path_sequence.append(id)
    out_edges = graph.out_edges(id)

    if len(out_edges) == 0:
        yield path_sequence
    else:
        for edge in out_edges:
            yield from path_recursion(graph=graph, alleles=alleles, id=edge[1], path_sequence=copy(path_sequence))


def enumerate_paths(alleles, graph):
    # Get start node
    start_id = None

    for i in range(len(alleles)):
        if alleles[i].is_left_flank:
            start_id = i

    print("Starting path recursion from %d" % start_id)

    paths = [p for p in path_recursion(graph=graph, alleles=alleles, id=start_id)]

    return paths


def test():
    path_weights = {
        (1,2,3):0,
        (4,5,6):1,
        (7,8,9):2,
        tuple([10]):3
    }

    paths = Paths()

    for p,[path,w] in enumerate(path_weights.items()):
        paths.increment_weight(path,w)
        paths.increment_weight(path,0)

        n2 = paths.get_path_name(path)
        p2 = paths.get_path_id(path)
        w2 = paths.get_path_weight(path)

        expected_name = '_'.join(list(map(str,path)))
        if not n2 == expected_name:
            raise Exception("ERROR: %s != %s" % (n2, expected_name))

        if not p2 == p:
            raise Exception("ERROR: %d != %d" % (p2, p))

        if not w2 == w:
            raise Exception("ERROR: %d != %d" % (w2, w))

        paths.increment_weight(path,10)
        paths.increment_weight(p,10)
        paths.increment_weight(expected_name,10)

    weights = paths.id_to_weight
    expected_weights = [30, 31, 32, 33]
    if not weights == expected_weights:
        raise Exception("ERROR: path weights %s != %s" % (str(weights), str(expected_weights)))

    print("SUCCESS")


if __name__ == "__main__":
    test()

