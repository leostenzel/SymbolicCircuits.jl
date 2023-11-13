abstract type G end

abstract type UHG <: G end
abstract type SG <: G end
abstract type RG <: G end

abstract type Hamiltonian <: G end

struct Pauli <: Hamiltonian
    mapping::Dict
end

struct gX <: UHG end
struct gY <: UHG end
struct gZ <: UHG end
struct gH <: UHG end

struct gT <: SG end
struct gS <: SG end


struct rX <: RG
    theta::Vector{Any}
end

struct rY <: RG
    theta::Vector{Any}
end

struct rZ <: RG
    theta::Vector{Any}
end


struct A <: Hamiltonian end
struct Ad <: Hamiltonian end


"""dagger gate"""
struct gXd <: UHG end
struct gYd <: UHG end
struct gZd <: UHG end
struct gHd <: UHG end

struct gTd <: SG end
struct gSd <: SG end

struct One end
struct Positive end
struct Negative end

struct rXd <: RG
    theta::Vector{Any}
end

struct rYd <: RG
    theta::Vector{Any}
end

struct rZd <: RG
    theta::Vector{Any}
end


import Base.(*)
function (*)(a::T, b::T) where {T<:RG}
    theta = copy(a.theta)
    append!(theta, b.theta)
    return T(theta)
end

abstract type Q end

struct Loc <: Q
    index::Int64
end

struct cLoc <: Q
    index::Int64
end

abstract type D end
struct Normal <: D end
struct Dagger <: D end

struct Gate{T<:D}
    g::G
    loc::Vector{Q}
end

const UGate = Gate{Normal}
const DaggerGate = Gate{Dagger}

function gate2string(g::Gate)
    type_dict = Dict(gX => "X", gY => "Y", gZ => "Z",
        gH => "H", gT => "T", gS => "S",
        rX => "RX", rY => "RY", rZ => "RZ",)

    g_str = type_dict[typeof(g.g)]
    for loc in g.loc
        if isa(loc, Loc)
            g_str = g_str * " " * string(loc.index)
        elseif isa(loc, cLoc)
            g_str = g_str * " " * "c" * string(loc.index)
        end
    end
    return g_str
end


loc_indices(a::Gate) = [x.index for x in a.loc]

function is_loc_intersect(a::Gate, b::Gate)
    indices1 = loc_indices(a)
    indices2 = loc_indices(b)
    is_intersect(indices1, indices2)
end

function is_loc_identity(a::Gate, b::Gate)
    indices1 = loc_indices(a)
    indices2 = loc_indices(b)
    is_set_identity(indices1, indices2)
end

function is_loc_type_identity(a::Gate, b::Gate)
    check = true
    for loca in a.loc
        for locb in b.loc
            if loca.index == locb.index
                if isa(loca, Loc)
                    check = check && isa(locb, Loc)
                elseif isa(loca, cLoc)
                    check = check && isa(locb, cLoc)
                end
            end
        end
    end
    return check
end

is_gate_type_identity(a::Gate, b::Gate) = typeof(a.g) == typeof(b.g)

function is_cancel(a::Gate, b::Gate)
    if is_loc_identity(a, b) && is_loc_type_identity(a, b) && is_gate_type_identity(a, b) && isa(a.g, UHG) || is_AAdagger(a, b)
        return true
    end
    return false
end

function is_gate_type_inverse(a::Gate, b::Gate)
    if isa(a.g, Hamiltonian) || isa(b.g, Hamiltonian)
        return false
    else
        if (a isa Gate{Normal} && b isa Gate{Dagger}) || (b isa Gate{Normal} && a isa Gate{Dagger})
            return true
        end
    end

    return false
end

function is_AAdagger(a::Gate, b::Gate)
    if is_loc_identity(a, b) && is_loc_type_identity(a, b) && is_gate_type_inverse(a, b)
        if isa(a.g, Union{UHG,SG})
            return true
        elseif a.g isa RG && b.g isa RG
            if is_set_identity(a.g.theta, b.g.theta)
                return true
            end
        else
            error()
        end
    end

    return false
end


function is_no_control(a::Gate)
    for loc in a.loc
        if loc isa cLoc
            return false
        end
    end
    return true
end

is_single_qubit(a::Gate) = length(a.loc) == 1

function is_Z(a::Gate)
    check = isa(a.g, gZ)
    check && is_single_qubit(a)
end

function is_X(a::Gate)
    check = isa(a.g, gX)
    check && is_single_qubit(a)
end


function is_Y(a::Gate)
    check = isa(a.g, gY)
    check && is_single_qubit(a)
end

function is_H(a::Gate)
    check = isa(a.g, gH)
    check && is_single_qubit(a)
end

function is_T(a::Gate)
    check = isa(a.g, gT)
    check && is_single_qubit(a)
end

function is_S(a::Gate)
    check = isa(a.g, gS)
    check && is_single_qubit(a)
end

function get_CNOT_loc_index(a::Gate, ::Val{:loc})
    for loc in a.loc
        if loc isa Loc
            return loc.index
        end
    end
end

function get_CNOT_loc_index(a::Gate, ::Val{:cloc})
    for loc in a.loc
        if loc isa cLoc
            return loc.index
        end
    end
end

function is_CNOT(a::Gate)
    if isa(a.g, gX)
        if length(a.loc) == 2
            if isa(a.loc[1], Loc) || isa(a.loc[2], Loc)
                if isa(a.loc[1], cLoc) || isa(a.loc[2], cLoc)
                    return true
                end
            end
        end
    end
    return false
end


function is_CNOT_T_commute(a::Gate, b::Gate)
    if (a isa Gate{Normal} && b isa Gate{Normal}) || (a isa Gate{Dagger} && b isa Gate{Dagger})
        if is_CNOT(a) && (is_T(b) || is_S(b) || is_Z(b))
            for loc in a.loc
                if isa(loc, cLoc)
                    index = loc.index
                    if index == b.loc[1].index
                        return true
                    end
                end
            end
        end
    end

    return false
end



function is_commute(a::Gate, b::Gate)
    @show a b
    if !is_loc_intersect(a, b)
        return true
    end

    if is_CNOT_T_commute(a, b) || is_CNOT_T_commute(b, a)
        return true
    end

    return false
end

function is_r_merge(a::Gate, b::Gate)
    if is_single_qubit(a) && is_single_qubit(b)
        if is_loc_identity(a, b) && is_gate_type_identity(a, b) && is_loc_type_identity(a, b)
            if isa(a.g, RG)
                return true
            end
        end
    end
    return false
end


function is_expand(a::Gate)
    if is_S(a)
        return true
    elseif is_Z(a)
        return true
    end
    return false
end



function expand(a::Gate)
    if is_S(a)
        index = a.loc[1].index
        g = typeof(a)(gT(), [Loc(index),])
        return :($(g) * $(g))

    elseif is_Z(a)
        index = a.loc[1].index
        g = typeof(a)(gS(), [Loc(index),])
        return :($(g) * $(g))

    end
    error()
end







function is_merge(a::Gate, b::Gate)
    @show a b
    if is_loc_identity(a, b)
        if is_S(a) && is_S(b) && !is_AAdagger(a, b)
            return true
        elseif is_T(a) && is_T(b) && !is_AAdagger(a, b)
            return true
        elseif is_r_merge(a, b)
            return true
        end
    end
    return false
end

function merge(a::Gate, b::Gate)
    if is_S(a)
        @show a.loc
        return :($(typeof(a))(gZ(), $(a.loc)))

    elseif is_T(a)
        return :($(typeof(a))(gS(), $(a.loc)))

    elseif is_r_merge(a, b)
        return :($(typeof(a))($(a.g * b.g), $(a.loc)))
    end
    error()
end


function to_cancel(a::Gate, b::Gate)
    if is_cancel(a, b)
        return :(One())
    else
        return :($a * $b)
    end
end


function to_dagger(a::Gate)
    if isa(a.g, Union{UHG,SG})
        if a isa Gate{Normal}
            g = Gate{Dagger}(typeof(a.g)(), a.loc)
        elseif a isa Gate{Dagger}
            g = Gate{Normal}(typeof(a.g)(), a.loc)
        else
            error()
        end

    elseif isa(a.g, RG)
        theta = a.g.theta
        if a isa Gate{Normal}
            g = Gate{Dagger}(typeof(a.g)(theta), a.loc)
        elseif a isa Gate{Dagger}
            g = Gate{Normal}(typeof(a.g)(theta), a.loc)
        else
            error()
        end
    else
        error()
    end

    return :($(g))
end


"""generate from vacuum"""
function generate_HXH(a::Gate)
    h = typeof(a)(gH(), a.loc)
    x = typeof(a)(gX(), a.loc)
    return :($(h) * $(x) * $(h))
end

function generate_HZH(a::Gate)
    h = typeof(a)(gH(), a.loc)
    z = typeof(a)(gZ(), a.loc)
    return :($(h) * $(z) * $(h))
end

function generate_HYH(a::Gate)
    h = typeof(a)(gH(), a.loc)
    y = typeof(a)(gY(), a.loc)
    return :($(h) * $(y) * $(h))
end

function is_HXH(a::Gate, b::Gate, c::Gate)
    return is_H(a) && is_X(b) && is_H(c) && is_loc_identity(a, b) && is_loc_identity(a, c) &&
           !is_gate_type_inverse(a, b) && !is_gate_type_inverse(a, c)
end

function is_HZH(a::Gate, b::Gate, c::Gate)
    return is_H(a) && is_Z(b) && is_H(c) && is_loc_identity(a, b) && is_loc_identity(a, c) &&
           !is_gate_type_inverse(a, b) && !is_gate_type_inverse(a, c)
end
function is_HYH(a::Gate, b::Gate, c::Gate)
    return is_H(a) && is_Y(b) && is_H(c) && is_loc_identity(a, b) && is_loc_identity(a, c) &&
           !is_gate_type_inverse(a, b) && !is_gate_type_inverse(a, c)
end

function is_XYZ(a::Gate, b::Gate, c::Gate)
    return is_X(a) && is_Y(b) && is_Z(c) && is_loc_identity(a, b) && is_loc_identity(a, c) &&
           !is_gate_type_inverse(a, b) && !is_gate_type_inverse(a, c)
end


function generate_Z(a::Gate)
    z = typeof(a)(gZ(), a.loc)
    return :($z)
end


function generate_X(a::Gate)
    x = typeof(a)(gX(), a.loc)
    return :($x)
end

function generate_Y(a::Gate)
    x = typeof(a)(gY(), a.loc)
    return :($x)
end

"""generate from vacuum end"""
