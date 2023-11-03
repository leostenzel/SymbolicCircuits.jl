using Metatheory
using SymbolicCircuits
using Test, TestItems


@testitem "test block merge" begin
    using SymbolicCircuits: is_block_merge, block_merge

    x1 = UGate(gX(), [Loc(1),])
    h2 = UGate(gH(), [Loc(2),])

    x23 = UGate(gX(), [Loc(2), Loc(3)])

    cnot_1c3 = UGate(gX(), [Loc(1), cLoc(3)])

    rx1 = UGate(rX([:theta1,]), [Loc(3),])


    b = block_merge(cnot_1c3, x1)
    # @show b
    @test !is_block_merge(x1, h2)
    @test is_block_merge(x23, h2)
    @test !is_block_merge(b, h2)
    @test !is_block_merge(b, rx1)
    # @show block_merge(b, rx1)
    # @show block_merge(rx1, b)
end


@testitem "test blocks" begin
    using Yao
    using SymbolicCircuits: merge_block2yao

    function test(block, circ, n_qubits)
        yao_gate = merge_block2yao(block, n_qubits)
        @show yao_gate

        yao_circ = to_yao(circ)

        a = zero_state(n_qubits) |> yao_gate
        b = zero_state(n_qubits) |> yao_circ

        @test a.state == b.state
    end

    x1 = UGate(gX(), [Loc(1),])
    x2 = UGate(gX(), [Loc(2),])
    z1 = UGate(gZ(), [Loc(1),])
    z2 = UGate(gZ(), [Loc(2),])
    y1 = UGate(gY(), [Loc(1),])
    h2 = UGate(gH(), [Loc(2),])
    y3 = UGate(gY(), [Loc(3),])
    cnot_1c3 = UGate(gX(), [Loc(1), cLoc(3)])

    block = ClauseBlock([x1, z2], [1, 2])
    circ = x1 * z2
    test(block, circ, 2)

    block = ClauseBlock([x1, h2], [1, 2])
    circ = x1 * h2
    test(block, circ, 2)

    block = ClauseBlock([y1, x2, z1, h2], [1, 2])
    circ = y1 * x2 * z1 * h2
    test(block, circ, 2)

    block = ClauseBlock([x1, y3, cnot_1c3], [1, 3])
    circ = x1 * y3 * cnot_1c3
    test(block, circ, 3)
end