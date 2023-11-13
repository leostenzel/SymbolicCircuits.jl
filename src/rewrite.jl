using Metatheory: Prewalk, Postwalk, PassThrough, Chain, Fixpoint

function dagger_rewriter()
    r = to_dagger_rule
    r = Postwalk(PassThrough(r))
    return r
end

function z2hxh_rewriter()
    r = Z2HXH_rule
    z2hxh = Postwalk(PassThrough(r))
    return z2hxh
end

function x2hzh_rewriter()
    r = X2HZH_rule
    z2hxh = Postwalk(PassThrough(r))
    return z2hxh
end


function block_simplify_rewriter()
    r = @rule a a::Block => block_simplify2expr(a)
    return Postwalk(PassThrough(r))
end