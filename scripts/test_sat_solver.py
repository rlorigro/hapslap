from ortools.sat.python import cp_model


def main():
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

    model.Add(b1 + b2 + b3 <= 2)

    # Objective
    # A case constructed so that every read has a different preferred haplotype.
    # The constraints on ploidy will need to prevent reads from being distributed evenly.
    model.Minimize(
        # r1
        3*r11a + 3*r11b +
        2*r12a + 2*r12b +
        1*r13a + 1*r13b +
        # r2
        1*r21a + 1*r21b +
        3*r22a + 3*r22b +
        2*r23a + 2*r23b +
        # r3
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
    main()
