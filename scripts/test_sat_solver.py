from ortools.sat.python import cp_model
from collections import defaultdict
import sys


def test2():
    samples = [1,2]
    reads = [1,2,3,4,5,6]
    paths = [1,2,3]
    path_to_read_costs = dict()
    path_to_read_vars = dict()
    path_to_sample_vars = dict()
    sample_to_read = {
        1:[1,2,3],
        2:[4,5,6]
    }

    # Sample A
    path_to_read_costs[(1,1)] = 3
    path_to_read_costs[(1,2)] = 2
    path_to_read_costs[(1,3)] = 1

    path_to_read_costs[(2,1)] = 1
    path_to_read_costs[(2,2)] = 3
    path_to_read_costs[(2,3)] = 2

    path_to_read_costs[(3,1)] = 1
    path_to_read_costs[(3,2)] = 2
    path_to_read_costs[(3,3)] = 3

    # Sample B
    path_to_read_costs[(3,4)] = 3
    path_to_read_costs[(1,4)] = 2
    path_to_read_costs[(2,4)] = 1

    path_to_read_costs[(3,5)] = 1
    path_to_read_costs[(1,5)] = 3
    path_to_read_costs[(2,5)] = 2

    path_to_read_costs[(3,6)] = 1
    path_to_read_costs[(1,6)] = 2
    path_to_read_costs[(2,6)] = 3

    model = cp_model.CpModel()

    # Define read assignment variables
    for edge in path_to_read_costs.keys():
        # 'edge' is a tuple with path_id,read_id
        path_to_read_vars[edge] = model.NewIntVar(0, 1, "p%dr%d" % edge)
        print("making variable: p%dr%d" % edge)

    # Constraint: each read must map to only one haplotype/path
    for read_id in reads:
        model.Add(sum([path_to_read_vars[(path_id,read_id)] for path_id in paths]) == 1)

    # Constraint: each sample's reads must map to at most two haplotypes/paths
    # Use a boolean indicator to tell whether any of a sample's reads are assigned to each haplotype/path
    for sample_id,read_group in sample_to_read.items():
        for path_id in paths:
            edge = (path_id,sample_id)
            path_to_sample_vars[edge] = model.NewBoolVar("p%ds%d" % edge)

            s = sum([path_to_read_vars[(path_id,read_id)] for read_id in read_group])
            model.Add(s >= 1).OnlyEnforceIf(path_to_sample_vars[edge])
            model.Add(s == 0).OnlyEnforceIf(path_to_sample_vars[edge].Not())

    # Now that the boolean indicators have been defined, use them to add a constraint on ploidy per sample
    for sample_id,read_group in sample_to_read.items():
        model.Add(sum([path_to_sample_vars[(path_id,sample_id)] for path_id in paths]) <= 2)

    # Cost function
    model.Minimize(sum([c*path_to_read_vars[e] for e,c in path_to_read_costs.items()]))

    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    assert status == cp_model.OPTIMAL

    print()
    print("=====Stats:======")
    print(solver.SolutionInfo())
    print(solver.ResponseStats())

    print("---- COSTS ----")
    for sample_id,read_group in sample_to_read.items():
        print("SAMPLE %d" % sample_id)

        for read_id in read_group:
            sys.stdout.write("\tr%d" % read_id)
        sys.stdout.write("\n")

        for path_id in paths:
            sys.stdout.write("p%d\t" % path_id)

            for read_id in read_group:
                c = path_to_read_costs[(path_id,read_id)]
                sys.stdout.write("%d\t" % c)

            sys.stdout.write("\n")

    print()
    print("---- RESULTS ----")
    for sample_id,read_group in sample_to_read.items():
        print("SAMPLE %d" % sample_id)

        for read_id in read_group:
            sys.stdout.write("\tr%d" % read_id)
        sys.stdout.write("\tis_path_used\n")

        for path_id in paths:
            sys.stdout.write("p%d\t" % path_id)

            for read_id in read_group:
                v = solver.Value(path_to_read_vars[(path_id,read_id)])
                c = path_to_read_costs[(path_id,read_id)]
                sys.stdout.write("%d\t" % v)

            sys.stdout.write("%d" % solver.Value(path_to_sample_vars[(path_id,sample_id)]))
            sys.stdout.write("\n")

    return


def test():
    model = cp_model.CpModel()

    # Variables

    # Sample A:
    # 3 reads, 3 haplotypes
    r11a = model.NewIntVar(0, 1, 'r11a')
    r12a = model.NewIntVar(0, 1, 'r12a')
    r13a = model.NewIntVar(0, 1, 'r13a')
    r21a = model.NewIntVar(0, 1, 'r21a')
    r22a = model.NewIntVar(0, 1, 'r22a')
    r23a = model.NewIntVar(0, 1, 'r23a')
    r31a = model.NewIntVar(0, 1, 'r31a')
    r32a = model.NewIntVar(0, 1, 'r32a')
    r33a = model.NewIntVar(0, 1, 'r33a')

    # Sample B:
    # 3 reads, 3 haplotypes
    r11b = model.NewIntVar(0, 1, 'r11b')
    r12b = model.NewIntVar(0, 1, 'r12b')
    r13b = model.NewIntVar(0, 1, 'r13b')
    r21b = model.NewIntVar(0, 1, 'r21b')
    r22b = model.NewIntVar(0, 1, 'r22b')
    r23b = model.NewIntVar(0, 1, 'r23b')
    r31b = model.NewIntVar(0, 1, 'r31b')
    r32b = model.NewIntVar(0, 1, 'r32b')
    r33b = model.NewIntVar(0, 1, 'r33b')

    # A read can only be assigned to one haplotype
    model.Add(r11a + r12a + r13a == 1)
    model.Add(r21a + r22a + r23a == 1)
    model.Add(r31a + r32a + r33a == 1)
    model.Add(r11b + r12b + r13b == 1)
    model.Add(r21b + r22b + r23b == 1)
    model.Add(r31b + r32b + r33b == 1)

    # Sample A ploidy constraints
    a1 = model.NewBoolVar("a1")

    model.Add(r11a + r21a + r31a >= 1).OnlyEnforceIf(a1)
    model.Add(r11a + r21a + r31a == 0).OnlyEnforceIf(a1.Not())

    a2 = model.NewBoolVar("a2")

    model.Add(r12a + r22a + r32a >= 1).OnlyEnforceIf(a2)
    model.Add(r12a + r22a + r32a == 0).OnlyEnforceIf(a2.Not())

    a3 = model.NewBoolVar("a3")

    model.Add(r13a + r23a + r33a >= 1).OnlyEnforceIf(a3)
    model.Add(r13a + r23a + r33a == 0).OnlyEnforceIf(a3.Not())

    model.Add(a1 + a2 + a3 <= 2)

    # Sample B ploidy constraints
    b1 = model.NewBoolVar("b1")

    model.Add(r11b + r21b + r31b >= 1).OnlyEnforceIf(b1)
    model.Add(r11b + r21b + r31b == 0).OnlyEnforceIf(b1.Not())

    b2 = model.NewBoolVar("b2")

    model.Add(r12b + r22b + r32b >= 1).OnlyEnforceIf(b2)
    model.Add(r12b + r22b + r32b == 0).OnlyEnforceIf(b2.Not())

    b3 = model.NewBoolVar("b3")

    model.Add(r13b + r23b + r33b >= 1).OnlyEnforceIf(b3)
    model.Add(r13b + r23b + r33b == 0).OnlyEnforceIf(b3.Not())

    # model.Add(b1 + b2 + b3 <= 2)

    # Objective
    # A case constructed so that every read has a different preferred haplotype.
    # The constraints on ploidy will need to prevent reads from being distributed evenly.
    model.Minimize(
        3*r11a + 3*r11b +
        2*r12a + 2*r12b +
        1*r13a + 1*r13b +
        #
        1*r21a + 1*r21b +
        3*r22a + 3*r22b +
        2*r23a + 2*r23b +
        #
        1*r31a + 1*r31b +
        2*r32a + 2*r32b +
        3*r33a + 3*r33b
    )

    # Solve
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    assert status == cp_model.OPTIMAL

    print("SAMPLE A")
    print("    %s %s %s %s" % ("r1", "r2", "r3", "h_active"))
    print("h1  %d  %d  %d  %d" % (solver.Value(r11a), solver.Value(r21a), solver.Value(r31a), solver.Value(a1)))
    print("h2  %d  %d  %d  %d" % (solver.Value(r12a), solver.Value(r22a), solver.Value(r32a), solver.Value(a2)))
    print("h3  %d  %d  %d  %d" % (solver.Value(r13a), solver.Value(r23a), solver.Value(r33a), solver.Value(a3)))

    print()
    print("SAMPLE B")
    print("    %s %s %s %s" % ("r1", "r2", "r3", "h_active"))
    print("h1  %d  %d  %d  %d" % (solver.Value(r11b), solver.Value(r21b), solver.Value(r31b), solver.Value(b1)))
    print("h2  %d  %d  %d  %d" % (solver.Value(r12b), solver.Value(r22b), solver.Value(r32b), solver.Value(b2)))
    print("h3  %d  %d  %d  %d" % (solver.Value(r13b), solver.Value(r23b), solver.Value(r33b), solver.Value(b3)))

    print()
    print("=====Stats:======")
    print(solver.SolutionInfo())
    print(solver.ResponseStats())


if __name__ == "__main__":
    test()
    # test2()
