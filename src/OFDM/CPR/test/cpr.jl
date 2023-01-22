
function cpr!(eq, r, ici_mtx, isi_mtx, X0, X1, flag0, flag1, rind, dind, D=size(eq,1); nvar=nothing, N0=nothing)
    Nt = size(X1,2); Nr = size(eq,2)
    maxk = maximum(dind) # データインデックスの最大値
    mink = minimum(dind) # データインデックスの最小値
    for q in 1:Nr
        for k in rind
            val = zero(ComplexF64) # 干渉係数
            ipow = zero(Float64) # 干渉電力
            for p in 1:Nt
                minr = k-D < mink ? mink : k-D
                maxr = k+D > maxk ? maxk : k+D
                for m in minr:maxr
                    c = ici_mtx[k,m,q,p]
                    s = isi_mtx[k,m,q,p]
                    if flag1
                        val += -c * X1[m,p]
                    end
                    if flag0
                        val += s * X0[m,p]
                    end
                    if !isnothing(nvar) && !isnothing(X1) && (m!==k)
                        ipow += real(abs2(c) * (1 - abs2(X1[m,p]))) # ICI残差電力
                    end
                    if !isnothing(nvar) && !isnothing(X0)
                        ipow += real(abs2(s) * (1 - abs2(X0[m,p]))) # ISI残差電力
                    end
                end
            end
            eq[k,q] = r[k,q] - val
            if !isnothing(nvar)
                nvar[k,q] = N0 + ipow
            end
        end
    end
end
