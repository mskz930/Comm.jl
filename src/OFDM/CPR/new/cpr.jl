using ...Digimod: SoftSymbol

# Cyclic-Property-Reconstruction
function cpr!(R, Y, Hisi, Hici, X1::T, X0::T, idxs; D=size(R,1)) where T<:Union{Nothing, AbstractArray{<:SoftSymbol}}
    nfft = size(R,1)
    for q = axes(R,2) # n_rx
        for k = inds # frequency
            val  = zero(eltype(Y)) # 干渉成分
            R[k,q] = Y[k,q]
            fst = ifelse(k-D<1, 1, k-D)
            lst = ifelse(k+D>nfft, nfft, k+D)
            if !isnothing(X0) # ISI
                for p = axes(X0,2) # n_tx
                    for m = fst:lst
                        val += Hisi[k,m,q,p] * X0[m,p].mean
                    end
                end
            end
            if !isnothing(X1) # ICI
                for p = axes(X1,2)
                    for m in fst:lst
                        if inds1[k,p] !== 0
                            val += -Hici[k,m,q,p]* X1[m,p].mean
                        end
                    end
                end
            end
            end
            R[k,q] -= val
        end
    end
end

#=
function cpr!(R, Y, Hisi, Hici, X1, X0, nvar, N0, inds, D)
    isi_flag, ici_flag = !isnothing(X1), !isnothing(X0)
    !isi_flag && !ici_flag && return
    nfft = size(R,1)
    for q = 1:n_rx
        for k = inds # subject index
            val  = zero(ComplexF64) # 干渉成分
            ipow = zero(Float64)    # 干渉電力
            R[k,q], nvar[k,q] = Y[k,q], N0
            fst = ifelse(k-D<1, 1, k-D)
            lst = ifelse(k+D>nfft, nfft, k+D)

            if isi_flag
                val, ipow += _isi_calc(Hisi, X0, fst, lst, q, k)
            end
            if ici_flag
                val, ipow += _ici_calc(Hici, X1, fst, lst, q, k)
            end
            R[k,q] -= val
            nvar[k,q] += ipow
        end
    end
end
=#


# ISI成分計算
function _isi_calc(Hisi, X0, fst, lst, k, q)
    val  = zero(eltype(Hisi))
    ipow = 0.0
    for p = axes(X0,2)
        for m = fst:lst
            coef = Hisi[k,m,q,p]
            val += coef * X0[m,p].mean
            ipow += real(abs2(coef) * X0[m,p].var) # ISI残差電力
        end
    end
    val, ipow
end

# ICI成分計算
function _ici_calc(Hici, X1, fst, lst, k, q)
    val  = zero(eltype(Hisi))
    ipow = 0.0
    for p = axes(X1,2)
        for m = fst:lst
            m == k && continue
            coef = Hici[k,m,q,p] # 干渉係数
            val += -coef * X1[m,p].mean
            (m!==k) && (ipow += real(abs2(c) * X1[m,p].var)) # ICI残差電力
        end
    end
    val, ipow
end
