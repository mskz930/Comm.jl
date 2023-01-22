include(joinpath(@__DIR__,"../MPA.jl"))

using Random
m,n = (10, 20)

# 適当な検査行列を作る
H = zeros(UInt8, m, n)
v = H[:,1]
v[1:3] .= 1
for n in axes(H,2)
    H[:,n] .= shuffle(v)
end

parity = MPA.Parity(H)
bnlist, cnlist  = MPA.makelist(parity.mat)

decoder = MPA.Decoder(parity, max_iter=1)

Lin = 10 .* randn(20)
Lout = decoder(Lin)

