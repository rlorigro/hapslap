from typing import Optional


class IncrementalIdMap:
    def __init__(self, offset=0):
        self.name_to_id = dict()
        self.id_to_name = [None for i in range(offset)]
        self.offset = offset
        self.n = offset

    def add(self, name: str) -> int:
        id = None

        if name not in self.name_to_id:
            id = self.n
            self.id_to_name.append(name)
            self.name_to_id[name] = id
            self.n += 1
        else:
            # Don't duplicate values if it already exists
            id = self.name_to_id[name]

        return id

    def get_id(self, name: str) -> Optional[int]:
        id = None

        if name in self.name_to_id:
            id = self.name_to_id[name]

        return id

    def get_name(self, id: int) -> Optional[str]:
        if not type(id) == int:
            raise Exception("ERROR: cannot get name for non-int id: " + str(id) + " of type " + str(type(id)))

        name = None

        if id < len(self.id_to_name):
            name = self.id_to_name[id]

        return name

    def write_to_file(self, output_path):
        with open(output_path, 'w') as file:
            file.write("id,name\n")
            for id,name in self:
                file.write("%s,%s\n"%(id,name))

    def get_count(self):
        return self.n - self.offset

    def get_max_id(self):
        return self.n - 1

    def ids(self):
        for id in range(self.offset,self.n):
            yield id

    def __len__(self):
        return self.get_count()

    def __iter__(self):
        for id in range(self.offset,self.n):
            yield id, self.id_to_name[id]


def test():
    items = [
        ["a",0],
        ["b",1],
        ["c",2]
    ]

    id_map = IncrementalIdMap()

    for name,id in items:
        id_map.add(name)

        id2 = id_map.get_id(name)
        if not id2 == id:
            raise Exception("ERROR: %d != %d" % (id2, id))

        name2 = id_map.get_name(id)
        if not name2 == name:
            raise Exception("ERROR: %s != %s" % (name2, name))

    items2 = list(id_map.name_to_id.items())
    for i in range(len(items)):
        if not tuple(items[i]) == tuple(items2[i]):
            raise Exception("ERROR: %s != %s" % (str(items2[i]), str(items[i])))

    print("SUCCESS")


def test_offset():
    offset = 100

    items = [
        ["a",100],
        ["b",101],
        ["c",102]
    ]

    id_map = IncrementalIdMap(offset=offset)

    for name,id in items:
        id_map.add(name)
        id_map.add(name)

        id2 = id_map.get_id(name)
        if not id2 == id:
            raise Exception("ERROR: %d != %d" % (id2, id))

        name2 = id_map.get_name(id)
        if not name2 == name:
            raise Exception("ERROR: %s != %s" % (name2, name))

    items2 = list(id_map.name_to_id.items())
    for i in range(len(items)):
        if not tuple(items[i]) == tuple(items2[i]):
            raise Exception("ERROR: %s != %s" % (str(items2[i]), str(items[i])))

    print("SUCCESS")


if __name__ == "__main__":
    test()
    test_offset()
