# crc.jl

"""
    CRC(inds)

CRC(巡回冗長検査ビット)
    
CRCの符号生成多項式の例)
CRC-4  : x^4 + x + 1

CRC-8  : x^8 + x^7 + x^6 + x^4 + x^2 + 1

CRC-16 : x^16 + x^12 + x^2 + 1

source: https://en.wikipedia.org/wiki/Cyclic_redundancy_check

"""
struct CRC{N}
    g::BitArray
end
# コンストラクタ
function CRC(inds)
    N = maximum(inds)
    g = BitArray(zeros(Bool, N+1))
    g[inds .+ 1] .= true
    CRC{N}(g)
end

# CRCのパリティビット計算
function (crc::CRC{N})(inp) where N
    g = crc.g
    res = BitArray(undef, N+1)
    inp = [inp; BitArray(zeros(N))]
    m = findfirst(isone, inp)-1 
    # @show m
    res[2:end] .= inp[m+1:m+N]
    while m < length(inp)-N
        circshift!(res, -1)
        res[end] = inp[m+N+1]
        # @show res
        # 先頭が0なら1になるまでシフトを繰り返す
        !res[1] && (m+=1; continue)
        for n in reverse(eachindex(g))
            res[end-n+1] ⊻= g[n]
        end
        m += 1
        # @show m, res
    end
    deleteat!(res, 1)
    res # parity bits
end