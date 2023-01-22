include("sorts.jl")

struct QRM_MLD{T}
end

function qrm_mld(y, H, refs, K)
    x = zeros(eltype(refs), size(R,2))
    qrm_mld!(x, y, H, refs, K)
end

function qrm_mld(z, R, refs, K)
    x = zeros(eltype(refs), size(R,2))
    qrm_mld!(x, z, Q, refs, K)
end

function qrm_mld!(x::AbstractVector{T}, y::AbstractVector{T}, H::AbstractMatrix{T}, refs, K) where T
    nrx, ntx = size(H)

    # SIR順に並び替えてQR分解する
    ord = columnsort(H, :SIR, rev=false)
    Hp = @view H[:,ord]
    Q, R = qr(Hp)

    qrm_mld!(x, z, R, order, refs, K)
end

function qrm_mld!()

    _m_algorithm()
end

# Mアルゴリズム
function _m_algorithm!(list, metrics, z, R, refs, M)
    T = eltype(refs)
    nref, ntx = length(refs), size(R,2)
    list_tmp    = Matrix{Int16}(undef, ntx, nref*M)
    metrics_tmp = Vector{Float64}(undef, nref*M)
    npath = 0
    # 探索開始
    l = 0
    while l < ntx
        if l == 0
            for i in eachindex(refs)
                metrics_tmp[i] = abs2(z[end] - (R[ntx,ntx] * refs[i]))
                list_tmp[ntx,i] = i
            end
            npath = ifelse(nref<M, nref, M)
            perm  = @views partialsortperm(metrics_tmp[1:nref], 1:npath)
            _save!(list, list_tmp, metrics, metrics_tmp, perm)
        else
            cnt = 0; m = ntx-l
            for i in 1:npath
                tmp = z[m]
                for j in m+1:ntx
                    tmp -= R[m,j] * refs[list[j,i]]
                end
                for j in 1:nref
                    cnt += 1
                    metrics_tmp[cnt] = metrics[i] + abs2(tmp - (R[m,m] * refs[j]))
                    list_tmp[m, cnt] = j
                    for k in m+1:ntx
                        list_tmp[k, cnt] = list[k, i]
                    end
                end
            end
            npath = ifelse((npath*nref)<M, npath*nref, M) # 生き残りパスの数
            perm = @views partialsortperm(metrics_tmp[1:cnt], 1:npath)
            _save_path!(list, list_tmp, metrics, metrics_tmp, perm)
        end
        l += 1
    end
    return list
end

# 生き残りパスの情報を記録する
function _save_path!(list, list_tmp, metrics, metrics_tmp, perm)
    for j in eachindex(perm)
        metrics[j] = metrics_tmp[perm[j]]
        for i in axes(list,1)
            list[i,j] = list_tmp[i,perm[j]]
        end
    end
end

