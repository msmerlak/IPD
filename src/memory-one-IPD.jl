using KrylovKit:eigsolve

## Press-Dyson determinant

# p = P(C|CC, CD, DC, DD)

function D(p, q, f)
    return @fastmath det([[-1 + p[1]*q[1], - 1 + p[1], - 1 + q[1], f[1]] [p[2]*q[3], -1 + p[2], q[3], f[2]] [p[3]*q[2], p[3], -1+q[2], f[3]] [p[4]*q[4], p[4], q[4], f[4]]])
end

function M(p, q)
    qq = [q[1], q[3], q[2], q[4]]
    return @. [p*qq p*(1-qq) (1-p)*qq (1-p)*(1-qq)]
end

function v(p, q)
    val, vec, conv  = eigsolve(transpose(M(p, q)), 1, :LR)
    return real(vec[1]./sum(vec[1]))
end