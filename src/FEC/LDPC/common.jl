
# ビット反転
function bitflip!(x)
    x .= -x
    return
end

# 符号
function mysign(x::Real)
    y = sign(x)
    if y == 0.0
        return 1.0
    else
        return y
    end
end

# atanh function
function myatanh(x::T) where T<:Number
    if -1.0 < x < 1.0
        return atanh(x)
    else
        return sign(x) * 19.07
    end
end

# phi function
function ϕ(x::T) where T <: Number
    y = abs(x)
    y = ifelse(y>38.1, 38.1, y)
    y = -log(tanh(y/2.0))
    return ifelse(isinf(y), 700.0, y)
end

# box-sum(soft-XOR)
function _boxsum(x,y)
    sign(x)*sign(y)*phi(phi(abs(x))+phi(abs(y)))
end


# 最大事後尤度比(LAPPs)の算出
function _lapp_calc!(Lapps::AbstractVector, Lch, rows, cnlist)
    for n in eachindex(Lch)
        Lapps[n] = Lch[n]
        neighbor = cnlist[n] # 周辺のチェックノード
        for i in neighbor
            Lapps[n] += rows[i][n]
        end
    end
    return Lapps
end

# 尤度比によるビット判定
function _decision!(dec_bits, Lapps)
    for i in eachindex(Lapps)
        dec_bits[i] = Lapps[i] < 0
    end
    return dec_bits
end

# チェックサムの計算
function checksum(Lpos, bnlist)
    cks = false
    for m in eachindex(bnlist) # Threads.@threads
        for j in bnlist[m]
            cks ⊻= (Lpos[j] < 0)
        end
        cks && break
    end
    cks
end
