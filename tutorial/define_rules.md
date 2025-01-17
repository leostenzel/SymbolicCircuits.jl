## How to define your own rules in `SymbolicCircuits.jl`?

One of the major goals of `SymbolicCircuits.jl` is to provide a symbolic system, such that users can define their own rules easily, and they can easily apply their rules to egraph provided by `Metatheory.jl` to show the effect of the rules on their circuit.

### What is the gate and circuit in `SymbolicCircuits.jl`?
To do this, firstly, let's take a look at the `Gate` struct provided by `SymbolicCircuits.jl`:

A gate struct is defined as below:
```julia
struct Gate
    g::G
    loc::Vector{Q}
end
```
Where `g` is an instance of type `G`, indicating which type of gate it is: Pauli gate? or rotate gate? or Hadamard?

`loc` is a Vector of instances of type `Q`, indicating its location on the circuit(which qubits it operate on?)

`abstract type G end` is an abstract type that has a bounch of subtypes:

    abstract type UHG <: G end where UHG stands for unitary & hermitian gates
        struct gX <:UHG end for Pauli X
        struct gY <:UHG end for Pauli Y
        struct gZ <:UHG end for Pauli Z
        struct gH <:UHG end for Hadamard
    abstract type SG <: G end where SG stands for other non-parametric(static) gates that are not UHG
        struct gS <:SG end for S gate
        struct gT <:SG end for T gate
    abstract type RG <: G end where RG stands for rotate gate
        struct rX <:RG theta::Vector{Any} end for rotate X gate where `theta` is just a vector of parameters it contains
        struct rY <:RG theta::Vector{Any} end for rotate Y gate
        struct rZ <:RG theta::Vector{Any} end for rotate Z gate

`SymbolicCircuits.jl` also defined dagger gate separately, for example, `struct gXd <:UHG end` stands for dagger X gate. However, currently I'm not sure if there are ways to avoid such a redundancy of definitions.

`abstract type Q end` is an abstract type which stands for qubit:

    struct Loc <: Q index::Int64 end for normal qubit the gate operates on
    struct cLoc <: Q index::Int64 end for control qubit

A single qubit pauli gate is define in the following way:
```julia
x1 = UGate(gX(), [Loc(1)])
```
And for CNOT operate on qubit 2 and controlled by qubit one:
```julia
x2c1 = UGate(gX(), [Loc(2), cLoc(1)])
```

rotate X with angle `theta1 + theta2` is defined:
```julia
rx1 = UGate(rX([:theta1, :theta2]), [Loc(1)])
```
The reason there allow multiple angles is in VQE system, a lot of ansatz has a linear mapping between the rotation angle and their parameters. A symbolic system like this will represent such phenomenon correctly.

As referred in the introductory document, the circuit in `SymbolicCircuits.jl` is defined as a chain of gate, connected by operator `*`, for example:
```julia
circ = x1 * x2c1 * rx1
```
It just looks like this, if you print it:
```
:(UGate(gX(), Q[Loc(1)]) * UGate(gX(), Q[Loc(2), cLoc(1)]) * UGate(rX(Any[:theta1, :theta2]), Q[Loc(1)]))
```

You see, it is just an expression!

The trick of `SymbolicCircuits.jl` is that it takes a chain of Ugate(quantum circuit) into an AST in julia, and systems like `Metatheory.jl` can easily handle it.

We are ready to define our own rules for quantum circuit :)

### the way to define the rules

To define a rule, we need to functions:

    a function to determine if the rule should be applied to an input object
    a function that manipulate the object

More specifically, consider a CNOT and T mutate rule:

![image](https://github.com/overshiki/SymbolicCircuits.jl/blob/main/tutorial/CNOT_T_commute.png)

We firstly defined a predicate function:
```julia
function is_CNOT(a::Gate)
    if isa(a.g, gX)
        if length(a.loc)==2
            if isa(a.loc[1], Loc) || isa(a.loc[2], Loc)
                if isa(a.loc[1], cLoc) || isa(a.loc[2], cLoc)
                    return true
                end
            end
        end
    end
    return false
end

function is_T(a::Gate)
    check = isa(a.g, gT)
    check = check && is_single_qubit(a)
    return check
end

function is_CNOT_T_commute(a::Gate, b::Gate)
    if is_CNOT(a) && is_T(b)
        for loc ∈ a.loc
            if loc isa cLoc
                if loc.index == b.loc[1].index
                    return true
                end
            end
        end
    end
    return false
end
```

And then generate the rule like this:
```julia
CNOT_T_commute_rule = @rule a b a::Gate * b::Gate => :($(b) * $(a)) where is_CNOT_T_commute(a, b)
```

That's it!

Note that the rule above just considered the case where `a` is CNOT and `b` is T, not vice verse way. Also, it does not consider the cases that, since CNOT and T has such a commute rule, so does CNOT and S, CNOT and Z. To do this, more complex rules related to expand S/Z into T and merge T into S/Z. For more details, one can refere to `src/rule.jl` on how to define them.
