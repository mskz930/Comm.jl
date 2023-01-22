module SphereDetector


using LinearAlgebra: qr, norm

#========== sphere decoding ============#
# include("./sphere_decoding.jl")

"""
    sphere_decoding(y, H, modparams, radius_init, dr=2) where T

MLdetection by sphere decoding algorithm
    - PAM/QAM modulated symbol is assummed.
"""
function sphere_decoding(y::AbstractVector{T}, Q::AbstractMatrix{T}, R::AbstractMatrix{T}, refs, slicer; radius_step=2.0) where T <: Real
    n_rx, n_tx = size(R)
    z = @views Q[:,1:n_tx]'*y[1:n_tx]

    x_ls = slicer(R \ z)
    srad = norm(z - R * x_ls)^2

    if n_rx > n_tx
        srad -= norm(Q[:,n_tx+1:end]'*y[n_tx+1:end])^2
    end
    # ２乗半径の初期値
    x_ml = _sphere_decoding(z, R, x_ls, refs, srad, radius_step)
    x_ml
end

function sphere_decoding(y::AbstractVector{T}, H::AbstractMatrix{T}, alphabets, slicer; srad=2.0, radius_step=2.0) where T <: Number
    n_rx, n_tx = size(H)
    if T <: Complex
        y = [real(y); imag(y)]
        H = [real(H) -imag(H);
            imag(H) real(H)]
    end

    Q, R = qr(H)
    x_ml, metric = sphere_decoding(y, Q, R, alphabets, slicer; radius_step=radius_step)
    return x_ml, metric
end

"""
Sphereリスト復号:
    候補点が1つ以上見つかるまで，半径を拡大し探索を繰り返す
    一つ見つかったあとは，トレースバックしてすべての候補ベクトルを探索する．
"""
function sphere_list_decoding(y::AbstractVector{T}, H::AbstractMatrix{T}, alphabets; srad=2.0, radius_step=2.0) where T
    nrx, ntx = size(H)
    Q, R = qr(H)
    R = R[1:ntx,:]
    z = @views Q[:,1:ntx]'*y[1:ntx]
    if nrx > ntx
        srad -= norm(Q[:,ntx+1]'*y[ntx+1:end])^2
    end
    x_cands, metrics = _sphere_list_decoding(z, R, alphabets, srad, radius_step)
    return x_cands, metrics
end

# Sphere復号の内部関数
# x_cand    : 現在のML候補点
# x_hat     : 候補点ベクトル
# ml_metric : 最尤メトリック
# deapth-first(深さ優先)のアルゴリズム
function _sphere_decoding(z, R, x_ls, X, srad, ϵ)
    m, n = size(R)
    cand_lists = [Int[] for _ in 1:n_tx] # 候補点のリスト
    metric = Inf
    x_cand = zero(x_ls)

    stage = 1
    flag = false
    while stage < n_tx+2
        if stage == 0 # 半径の更新ｎ
            flag && break
            srad *= ϵ
            stage += 1
        elseif 1 <= stage <= n_tx # 現在の超球半径での探索
            stage = _ml_search(z, R, x_cand, x_ls, cand_lists, X, srad, stage)
        elseif stage == n+1 # 現在の超球内に候補点のみであることを確認
            metric, srad, stage = _is_ml_solution(z, R, x_cand, x_ls, cand_lists, metric, srad, stage)
            flag = true
        else
            error("dim=$dim<0!")
        end
    end
    x_hat, ml_metric
end


# 各ステージでの処理
function _ml_search(z, R, x_cand, x_ls, cand_lists, X, srad, layer)
    n = size(R,2)
    ϵ = 1e-3
    i = n - layer + 1 # 探索位置index

    # 境界を計算し，候補点を加える
    cand_list = cand_lists[i] # 現階層の候補点集合
    lower_bound, upper_bound = _calc_bounds(z, R, x_cand, x_ls, X, srad, layer)

    # 境界条件を満たすシンボル候補を保存する
    for (idx,x) in enumerate(X)
        if lower_bound-ϵ <= x <= upper_bound+ϵ
            push!(cand_list, idx)
        end
    end

    if !isempty(cand_list)
        # 境界を満たすシンボルが見つかった場合
        # => 候補シンボルの一番上を取り出して
        # 候補シンボルベクトルのi番目要素を代入する
        idx = popfirst!(cand_list)
        x_cand[i] = X[idx]
        layer += 1
    else
        # 境界を満たす候補点が見つからなかった場合
        layer -= 1
        while layer > 1
            i = n-layer+1
            if isempty(cand_lists[i])
                layer -= 1
            else
                idx = popfirst!(cand_lists[i])
                x_cand[i] = X[idx]
                layer += 1
            end
        end
    end
    return layer
end

# C0: squared radius
function _calc_bounds(z, R, x_cand, x_ls, Χ, srad, layer)
    n = size(R,2) # ベクトル長
    i = n - layer + 1 # anteanna index

    # sqrt term
    α_tmp = 0.0
    for k in n:-1:i+1
        tmp = 0.0
        for j in k:n
            tmp += R[k,j] * (x_cand[j] - x_ls[j])
        end
        α_tmp += abs2(tmp)
    end
    α = srad - α_tmp

    # no sqrt term
    β = 0.0
    for j in i+1:n # current layer
        β += R[i,j] * (x_cand[j] - x_ls[j])
    end

    upper_bound = x_ls[i] + (α - β) / abs(R[i,i])
    lower_bound = x_ls[i] - (α - β) / abs(R[i,i])

    lower_bound, upper_bound
end


# 候補点ベクトルの比較
function _is_ml_solution(z, R, x_ml, x_cand, x_ls, cand_lists, metric, srad, stage)
    # 候補点との二乗距離を計算
    metric = _squared_radius(z .- R * x_cand)

    # 候補点リストが空
    # => 現在の半径における超級内に存在する候補点ベクトルはx_cand以外存在しない
    if _isempty(cand_set)
        flag = 1
    end

    # 半径を変えて候補点を探す
    if isinf(metric)
        srad = ml_metric = _squared_radius(z, R, x_cand)
        _empty!(cand_set)

    end
    if ml_metric > x_metric
        x_ml .= x_cand
        metric = srad
        _empty!(cand_set)
    else
        ndim = length(x_cand)
        stage -= 1
        while stage > 1
            idx = ndim-stage+1
            if isempty(cand_lists[idx])
                stage -= 1
            else
                x_cand[idx] = popfirst!(cand_lists[idx])
                stage += 1
            end
        end
    end
    return metric, srad, stage, flag
end


# 候補点集合からリスト長を返す
function _list_size(list_set)
    if isempty(list_set)
        return 0
    else
        list_size = 1
        for elm in list_set
            list_size *= length(elm)
        end
        return list_size
    end
end

# 距離の二乗
function _squared_radius(z, R, x̂)
    r = z - R * x̂
    return dot(r,r)
end

# リストの更新
function _list_update!(list, col_id)
    @views list[1:end-1, col_id] .= list[2:end, col_id]
    list[end, col_id] = 0.0
    list
end

# リスト長の取得
function _get_list_size(list::Matrix{T}, col_id::Integer) where T
    list_size = 0
    for i in axes(list, 1)
        val = list[i, col_id]
        if val > 0
            list_size += 1
        else
            break
        end
    end
    list_size
end

# リストが唯一であるかを確認
function _isempty(cand_set)
    for elm in cand_set
        !isempty(elm) && return false
    end
    true
end

function _empty!(cand_set)
    for elm in cand_set
        !isempty(elm) && empty!(elm)
    end
end


# 前段の候補シンボルを変更する(なければ更に戻る)
function _trace_back(x_now, x_list, dim)
    while dim > 0
        list_size = _get_list_size(x_list, dim) # リストサイズの取得
        if list_size > 0
            x_now[stage] = x_list[1, stage]
            x_list = _list_update!(x_list, dim)
            dim += 1
            break
        else
            dim -= 1
        end
    end
    state = Int(dim == 0)
    state, dim
end



#=
# 境界計算: argmin x in |z - Rx|^2 <=  √C
function _bounds_calc(x_now, x_hat, R, srad, scale, stage)
    α = srad; β = 0.0
    ndim = length(x_hat)    # 候補点次元
    idx = ndim - stage + 1

    # alpha, beta 計算
    for i in idx+1:ndim
        temp_abs = 0.0
        for j in i:ndim
            temp = R[i,j] * (x_now[j] - x_hat[j])
            temp_abs += temp
        end
        α -= abs2(temp_abs)
    end
    α = sqrt(α)
    i = idx
    for j in i+1:ndim
        β += R[i,j] * (x_now[j]-x_hat[j])
    end

    lower_bound = x_hat[idx] - sign(R[idx,idx]) * (α + β)/R[idx,idx]
    upper_bound = x_hat[idx] + sign(R[idx,idx]) * (α - β)/R[idx,idx]
    lower_bound = lower_bound > 0 ? floor(lower_bound * scale) / scale : ceil(lower_bound * scale) / scale
    upper_bound = upper_bound > 0 ? floor(upper_bound * scale) / scale : ceil(upper_bound * scale) / scale

    return upper_bound, lower_bound
end
=#

# Sphereリスト復号の内部関数
function _sphere_list_decoding(z, R, alphabets, srad, radius_step)
    ndim = size(R,2)
    n_stage = ndim
    cand_list = [Int[] for _ in 1:ndim] # 候補点のリスト
    metrics = Float64[]
    x_cands    = Vector{Float64}[]
    x_hat     = Array{Float64}(undef, ndim)

    stage = 1
    flag = false
    while stage < ndim+2
        if stage == 0
            flag && break
            srad *= radius_step
            stage += 1
        elseif stage <= ndim # 現在の超球半径での探索
            stage = _stage_processing(x_hat, cand_list, z, R, srad, alphabets, stage)
        elseif stage == ndim+1 # 現在の超球内に候補点のみであることを確認
            srad, stage = _save_list!(metrics, x_cands, x_hat, cand_list, z, R, srad, stage)
            flag = true
        else
            error("dim=$dim<0!")
        end
    end
    x_cands, metrics
end

# 候補点ベクトルをリストに保存し，次の候補点の探索の準備をする
function _save_list!(metrics, x_cands, x_hat, cand_list, z, R, srad, stage)
    x_metric = _squared_radius(z, R, x_hat)

    push!(metrics, x_metric)
    push!(x_cands, copy(x_hat))

    ndim = length(x_hat)
    stage -= 1
    while stage > 0 # trace back
        idx = ndim-stage+1
        if isempty(cand_list[idx])
            stage -= 1
        else
            x_hat[idx] = popfirst!(cand_list[idx])
            return srad, stage+1
        end
    end
end

end # module
