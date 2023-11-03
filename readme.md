
# Symbolic system for Quantum computing

### Install

To install, type `]` in a julia (>=1.5) REPL and then input
```julia pkg
pkg> add https://github.com/leostenzel/SymbolicCircuits.jl
```

### Quick start
To define and simplify a quantum circuit, just do
```julia
using SymbolicCircuits
x1 = UGate(gX(), [Loc(1), ])
h2 = UGate(gH(), [Loc(2), ])
y3 = UGate(gY(), [Loc(3), ])
circ = x1 * x1 * h2 * y3 * x1 * x1 * h2
ncirc = egraph_simplify(circ, Val(:default_rule))
```
It will simplify `circ` into `ncirc` according to commute and cancel rules.

To test if two circuits are in equivalent according to commute and cancel rules, just do
```julia
SymbolicCircuits.areequal(Val(:default_rule), circ, ncirc)
```
It will return `true` if `circ` and `ncirc` are in equivalent, otherwise `missing`.

#### Working with Yao.jl
`SymbolicCircuits.jl` also provides method to turn symbolic circuit defined above to `YaoBlock` in `Yao.jl`, so that the full power of `Yao.jl` would be available
```
yao_circ = to_yao(circ; num_qubits=3)
println(yao_circ)
```
To plot the circuit, just do
```
to_yaoplot("./test/plot.png", circ)
```
It will plot the circuit using `YaoPlots.jl`

![image](./test/plot.png)


### What does `SymbolicCircuits.jl` provide?

`SymbolicCircuits.jl` provides a symbolic system for representation of Quantum circuit, in which, one can manipulate Quantum circuit using term rewriting & equality saturation techniques. Using this package, one can easily define any syntactic rules of Quantum circuit(for example, mutation between two quantum gates), and apply it to term rewriting and equality saturation Modules provided by[ `Metatheory.jl` ](https://github.com/JuliaSymbolics/Metatheory.jl)(Yes, This project is highly dependent on and highly motivated by [`Metatheory.jl`](https://github.com/JuliaSymbolics/Metatheory.jl)). Doing so, tasks such as circuit simplification, equivalence detection, code generation become easily achievable.

### What is the symbolic system of `SymbolicCircuits.jl`?
#### Quantum gates
In `SymbolicCircuits.jl`, quantum gates are represented by the instances of `Gate` struct, for example, the Pauli X gate on 1st Qubit is defined by
```julia
using SymbolicCircuits
x1 = UGate(gX(), [Loc(1), ])
```
Where `UGate` is alias for `Gate{Normal}`, where the counterpart is `Gate{Dagger}` which stands for `daggered` version of `UGate`.
`gX()` represents the type of Pauli X(gX is short for gate X), `Loc(1)` indicates the gate is applied on Qubit 1.
We can also define multi-Qubits gate, for example, CNOT gate on 1st Qubit and controlled by 2nd Qubit
```julia
cnot_1c2 = UGate(gX(), [Loc(1), cLoc(2)])
```
Where `cLoc(2)` indicates the gate is controlled by 2nd Qubit

We can also define parametric gate such as rotate X, controlled by parameter `:theta`
```julia
rx = UGate(rX([:theta, ]), [Loc(3), ])
```
Since simplifying many VQE circuits will often result in rotate gate controlled by multiple parameters(sum of multiple rotate angles), we also allow this
```julia
rx2 = UGate(rX([:theta1, :theta2, :theta3]), [Loc(3), ])
```
#### Quantum circuit
It is wellknown that, symbolically, quantum circuit is just a chain of quantum gate, applied from left to the right. In `SymbolicCircuits.jl`, we define quantum circuit as expression of gates connected by `*` operator.
```julia
using SymbolicCircuits
x1 = UGate(gX(), [Loc(1), ])
h2 = UGate(gH(), [Loc(2), ])
y3 = UGate(gY(), [Loc(3), ])
circ = x1 * x1 * h2 * y3 * x1 * x1 * h2
```
This allows flexible ways of defining a circuit, for example, we can also
```julia
using SymbolicCircuits
let circ = head_circuit()
    for _ in 1:2
        for g in [x1, h2, y3]
            circ *= g
        end
    end
end
```
#### Syntactic rules
The core idea of `SymbolicCircuits.jl` is that quantum circuit can be represented using symbolic expression and some rules(such as mutation rules) can be represented using syntactic rules of the expressions. Then term rewriting & equality saturation system could be used to manipulate the circuit. In `SymbolicCircuits.jl`, we provide easy ways to define different type of syntactic rules one may imagine. For example, a simple commute rule and cancel rule can be defined
```julia
using SymbolicCircuits
using Metatheory

#commute rules
com_rule = @rule a b a::Gate * b::Gate => :($(b) * $(a)) where is_commute(a, b)

#cancel rules
can_rule = @rule a b a::Gate * b::Gate => One() where is_cancel(a, b)

```
Where `is_commute` and `is_cancel` are functions provided in `src/gate.jl` to determine if two gates are commute and if they could be cancelled out.
Currently, `is_commute` considers two cases:
  - gate `a` and gate `b` do not have common `Loc` or `cLoc`
  - gate `a/b` is a `Z|S|T` gate, and gate `b/a` is a `CNOT` gate, where `cLoc` of `b/a` has the same index with `Loc` of `a/b`

`is_cancel` considers two cases:
  - gate `a` and gate `b` are identical and they belong to unitary & Hermitian gate
  - gate `a/b` is the dagger version of gate `b/a`

Note that the rule definition process here is just a normal rule definition process in [`Metatheory.jl`](https://github.com/JuliaSymbolics/Metatheory.jl)(users can also define rewrite rule using `-->` and equality rule using `==`, following the document of [`Metatheory.jl`](https://github.com/JuliaSymbolics/Metatheory.jl)), what `SymbolicCircuits.jl` provides is a system of gate expression and circuit expression where[ `Metatheory.jl` ](https://github.com/JuliaSymbolics/Metatheory.jl)could be applied to.

Additionally, `SymbolicCircuits.jl` also provides some built-in rules for applications including circuit simplification, equivalence detection(and code generation in the next release). They are included in `src/rule.jl`. Users of `SymbolicCircuits.jl` could of course define more rules they would like to. For more details, one can refer to the [tutorial](https://github.com/overshiki/SymbolicCircuits.jl/blob/main/tutorial/define_rules.md) on how to do this.

#### Term rewriting & Equality saturation
The core operation of `SymbolicCircuits.jl` is to apply a variety of syntactic rules to term rewriting & equality saturation system provided by [`Metatheory.jl`](https://github.com/JuliaSymbolics/Metatheory.jl), which allows powerful circuit manipulate functions. For example, term rewriting could be used to transform all Pauli Z gates in a circuit into `H X H` following rewriting rule `Z=HXH`. This could be easily handled using `SymbolicCircuits.jl`
```julia
using SymbolicCircuits
using SymbolicCircuits: is_Z, is_commute, is_cancel
using Metatheory
using Metatheory: PassThrough, Postwalk

x1 = UGate(gX(), [Loc(1), ])
y3 = UGate(gY(), [Loc(3), ])
z2 = UGate(gZ(), [Loc(2), ])
z3 = UGate(gZ(), [Loc(3), ])

circ = x1 * z2 * y3 * z3 * x1 * z3

function get_HXH(a::Gate)
    h = UGate(gH(), a.loc)
    x = UGate(gX(), a.loc)
    return :($(h) * $(x) * $(h))
end

r = @rule a a::Gate => get_HXH(a) where is_Z(a)
r = Postwalk(PassThrough(r))
@show Circuit(r(circ.expr))
```
Of course, this is just a simple rewriting rule, and could be done using any non-symbolic systems such as direct coding in qiskit|Yao.jl|mindquantum|...(for example, in mindquantum, circuits are represented as object of List class in python, such manipulation could be easily done using pop and insert operations).

The full power of `SymbolicCircuits.jl` comes from the equality saturation techniques in symbolic & compiler community. In equality saturation(Eqsat), the term rewriting are allowed to be handled without an order. For example, bidirectional rules such as gate commute rule `a*b==b*a` are allowed. The key point of Eqsat is that it allows all combination of rewriting operations to be stored in memory(usually in efficient EGraph data structure as implemented in [`Metatheory.jl`](https://github.com/JuliaSymbolics/Metatheory.jl)), and efficiently traversed. Thus, one could set the search goal to search for specific expressions in all equivalent expressions resulted from the rewriting rules, which would be a complex task for non-symbolic systems without Eqsat techniques, such as qiskit|Yao.jl|mindquantum|...

To demonstrate, let's consider circuit simplification operation. We consider only a combination of commute rule `a*b==b*a where a, b commute` and cancel rule `a*b->0 where a, b could be cancelled`. These simple rules already provide a strategy to simplify a circuit, i.e., for any pair of gates in a circuit, commute if possible, cancel if possible, and then try more commutation and cancellation, until no progress could be made. In a non-symbolic system without Eqsat, such strategy could only be implemented using a heuristic approach without a saturation guarantee(or it could have saturation guarantee but need a lot of codings). In `SymbolicCircuits.jl` it would be easily achievable and has a guarantee that all possible combination are traversed.
```julia
using SymbolicCircuits
using Metatheory
using Metatheory.Library: @right_associative, @left_associative
v = AbstractRule[]
push!(v, @rule a b a::Gate * b::Gate => :($(b) * $(a)) where is_commute(a, b))
push!(v, @rule a b a::Gate * b::Gate => One() where is_cancel(a, b))
push!(v, @rule a b b::One * a::Gate => :($(a)))
push!(v, @rule a b a::Gate * b::One => :($(a)))

ra = @right_associative (*)
la = @left_associative (*)
push!(v, ra)
push!(v, la)

function simplify(circuit)
    g = EGraph(circuit)
    params = SaturationParams(timeout=10, eclasslimit=40000)
    report = saturate!(g, v, params)
    circuit = extract!(g, astsize)
    return circuit
end

circ = x1 * x1 * h2 * y3 * x1 * x1 * h2
ncirc = simplify(circ)
```
result:
```julia
UGate(gY(), Q[Loc(3)])
```

Another powerful function provided by the system is that, it allows to detect if two circuit are in equivalence under certain rules. For example, continue from the above circuit simplification task, just do
```julia
areequal(v, circ, ncirc)
```
results in:
```julia
true
```
will tell if `circ` and `ncirc` are in equivalent under the rules `v`.


#### Practical usage
For those want to use `SymbolicCircuits.jl` as openbox toolset for some practical task, `SymbolicCircuits.jl` currently provides built-in functions including: circuit simplification, equivalence detection, where a default set of syntactic rules is used. (more simplification target and code generation function are on the way)
The list of the functions are:

`egraph_simplify`: circuit simplification using default set of rules
```julia
circ = x1 * x1 * h2 * y3 * x1 * x1 * h2
circ = egraph_simplify(circ, Val(:default_rule); verbose=false)
```
`areequal`: overload of `areequal` function in [`Metatheory.jl`](https://github.com/JuliaSymbolics/Metatheory.jl), using default set of rules
```julia
SymbolicCircuits.areequal(Val(:default_rule), circ1, circ2, circ3)
```

#### More information
For more information, please refer to `test` folder. Issue and PR are welcome.

#### Next step:
    - Some examples in example folder
    - More rewriting rules
    - T-gate reduction functions
    - Codegen for numerical simulation
    - Documentations!