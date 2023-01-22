module Polar

using Memoize

export PolarEncoder, PolarDecoder

""" Polar符合化オブジェクト
"""
struct Encoder
    N::Integer
end
PolarEncoder = Encoder


""" Polar復号化オブジェクト
"""
struct Decoder
end
PolarDecoder = Decoder



"""
    encode!(x, N=length(x))

Polar符合化: 入力バイナリxを再帰的に変換する

"""
function encode!(x; debug=false)
    N = length(x)
    N == 1 && return x

    nstep = 2
    while nstep <= N
        mstep = nstep >> 1 # XORの間隔
        # @show nstep, mstep

        k = mstep # 分割ブロック数
        for n in 1:k
            for m in n:nstep:N-1
                # @show m, m+mstep
                x[m] ⊻= x[m+mstep] # XOR
            end
        end
        debug && @show x
        nstep <<= 1 # nstep*2
    end
    x
end
encode(x) = encode!(copy(x))

"""
Polar復号
"""
function decode!()

end


end # module
