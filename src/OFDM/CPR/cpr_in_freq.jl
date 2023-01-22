

# CPR in Frequency domain
function cpr!(R, Y, Hisi, Hici, X1, X0; target_idxs=axes(R,1), D=div(size(R,1),2))
    nfft = size(R,1)
    for q = axes(R,2) # n_rx
        for k in target_idxs # freq
            left  = k-D:k-1
            right = k:k+D

            ival  = zero(eltype(Y)) # 干渉成分
            if !isnothing(X0) # ISI
                ival = _calc_isi_term(ival, Hisi, X0, left, right, k, q)
            end
            if !isnothing(X1) # ICI
                ival = _calc_ici_term(ival, Hici, X1, left, right, k, q)
            end
            R[k,q] = Y[k,q] - ival
        end
    end
    R
end


function cpr!(R, Y, Hisi, Hici, X1, flag1, X0, flag0, N, N0; target_idxs=axes(R,1), D=size(R,1))
    isi_flag, ici_flag = !isnothing(X1), !isnothing(X0)
    nfft = size(R,1)
    for q = axes(R,2)
        for k in target_idxs # subject index
            ival, ipow = zero(ComplexF64), zero(Float64)
            fst = ifelse(k-D<1, 1, k-D)
            lst = ifelse(k+D>nfft, nfft, k+D)
            if flag0
                ival, ipow = _isi_calc(ival, ipow, Hisi, X0, fst, lst, k, q)
            end
            if flag1
                ival, ipow = _ici_calc(ival, ipow, Hici, X1, fst, lst, k, q)
            end
            R[k,q] = Y[k,q] - ival
            N[k,q] = N0 + ipow
        end
    end
end


# ISI成分計算
function _calc_isi_term(ival, Hisi, X0::T, rng1, rng2, k, q) where T<:AbstractArray{<:Number}
    nfft = size(Hisi,1)
    for p = axes(X0,2)
        for m in flatten((rng1, rng2))
            m = ((m-1) + nfft) % nfft + 1
            ival += Hisi[k,m,q,p] * X0[m,p]
        end
    end
    ival
end
function _calc_isi_term(ival, Hisi, X0::T, rng1, rng2, k, q) where T<:AbstractArray{<:SoftSymbol}
    nfft = size(Hisi,1)
    for p = axes(X0,2)
        for m in flatten((rng1, rng2))
            m = ((m-1) + nfft) % nfft + 1
            ival += Hisi[k,m,q,p] * X0[m,p].mean
        end
    end
    ival
end
function _calc_isi_term(ival, ipow, Hisi, X0::T, rng1, rng2, k, q) where T<:AbstractArray{<:SoftSymbol}
    nfft = siz(Hisi,1)
    for p = axes(X0,2)
        for m in flatten((rng1,rng2))
            m = ((m-1) + nfft) % nfft + 1
            coef = Hisi[k,m,q,p]
            ival += coef * X0[m,p].mean
            ipow += real(abs2(coef) * X0[m,p].var) # ISI残差電力
        end
    end
    ival, ipow
end
function _calc_isi_term(ival, ipow, Hisi, X0::Nothing, fst, lst, k, q)
    nfft = size(Hisi,1)
    ival  = zero(eltype(Hisi))
    ipow = 0.0
    for p = axes(Hisi,4)
        for m in flatten((rng1,rng2))
            m = ((m-1) + nfft) % nfft + 1
            ipow += real(abs2(Hisi[k,m,q,p])) # ISI残差電力
        end
    end
    ival, ipow
end

# ICI成分計算
function _calc_ici_term(ival, Hici, X1::T, rng1, rng2, k, q) where T <: AbstractArray{<:Number}
    nfft = size(Hici,1)
    for p = axes(Hici,4)
        for m in flatten((rng1,rng2))
            m = ((m-1) + nfft) % nfft + 1
            m == k && continue
            ival += -Hici[k,m,q,p]* X1[m,p]
        end
    end
    ival
end
function _calc_ici_term(ival, Hici, X1::T, rng1, rng2, k, q) where T <: AbstractArray{<:SoftSymbol}
    nfft = size(Hici,1)
    for p = axes(Hici,4)
        for m in flatten((rng1,rng2))
            m = ((m-1) + nfft) % nfft + 1
            m == k && continue
            ival += -Hici[k,m,q,p]* X1[m,p].mean
        end
    end
    ival
end
function _calc_ici_term(ival, ipow, Hici, X1::T, rng1, rng2, k, q) where T <: AbstractArray{<:SoftSymbol}
    nfft = size(Hici,1)
    for p = axes(X1,2)
        for m in flatten((rng1,rng2))
            m = ((m-1) + nfft) % nfft + 1
            m == k && continue
            coef = Hici[k,m,q,p] # 干渉係数
            ival += -coef * X1[m,p].mean
            ipow += real(abs2(coef) * X1[m,p].var) # ICI残差電力
        end
    end
    ival, ipow
end
function _calc_ici_term(ival, ipow, Hici, X1::Nothing, rng1, rng2, k, q)
    nfft = size(Hici,1)
    for p = axes(Hici,4)
        for m in flatten((rng1,rng2))
            m = ((m-1) + nfft) % nfft + 1
            m == k && continue
            ipow += real(abs2(Hici[k,m,q,p])) # ICI残差電力
        end
    end
    ival, ipow
end


function cpr_in_freq!()

end
