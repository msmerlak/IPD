using DrWatson
@quickactivate

include("../_research/IPD.jl")

## Eq (5)
using KrylovKit

RSTP = [3., 0., 5., 1.]

f = rand(4)
p = q = RAND
D(p, q, f)


D(p, q, f)/D(p, q, ones(4)) ≈ dot(v(p, q), f)


p = q = TFT
D(p, q, f)/D(p, q, ones(4))
dot(v(p, q), f)


p = TFT
q = RAND
D(p, q, f)/D(p, q, ones(4)) ≈ dot(v(p, q), f)

p = RAND
q = TFT
D(p, q, f)/D(p, q, ones(4)) ≈ dot(v(p, q), f)

p = RAND
q = TFT
D(q, p, f)/D(q, p, ones(4)) ≈ dot(v(p, q), f)


f = rand(4)
D(RAND, RAND, f)

p = RAND
q = tft
D(p, q, f)/D(p, q, ones(4))
dot(v(p, q), f)


function minor(A, i, j)
    m, n = size(A)
    B = similar(A, m-1, n-1)
    for j′=1:j-1, i′=1:i-1; B[i′,j′]=A[i′,j′]; end
    for j′=1:j-1, i′=i+1:m; B[i′-1,j′]=A[i′,j′]; end
    for j′=j+1:n, i′=1:i-1; B[i′,j′-1]=A[i′,j′]; end
    for j′=j+1:n, i′=i+1:m; B[i′-1,j′-1]=A[i′,j′]; end
    return B
end

function adj(A)
    C = similar(A)
    n, m = size(A)
    for i = 1:n, j = 1:m
        C[i,j] = (-1)^(i+j)*det(minor(A, i, j))
    end
    return transpose(C)
end



