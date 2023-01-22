module SphereDetector

using LinearAlgebra: qr, dot


"""
    sphere_detector(y, H, modparams, radius_init, dr=2) where T

MLdetection by sphere decoding algorithm
    - PAM/QAM modulated symbol is assummed.
"""
function sphere_detector(y::AbstractVector{T}, H::AbstractMatrix{T}, refs, slicer; step_size=2.0, output=:ml) where T <: Number
    if T <: Complex
        y = [real(y); imag(y)]
        H = [real(H) -imag(H);
            imag(H) real(H)]
    end
    Q, R = qr(H)
    x_ml, metric, metrics, visited_list = sphere_detector(y, Q, R, refs, slicer, step_size)
    
    if output == :ml
        return x_ml, metric
    else # :list
        return @views metrics, visited_list
    end
end
function sphere_detector(y::AbstractVector{T}, Q::AbstractMatrix{T}, R::AbstractMatrix{T}, refs, slicer, step_size=2.0) where T <: AbstractFloat
    _, n = size(R)
    z = Q'*y
    x_zf = (R \ y)[1:n] |> slicer
    srad = _squared_error(z, R, x_zf) # 半径初期化
    x_ml = _sphere_detector(z, R, refs, srad, step_size)
    x_ml
end


# Sphere復号の内部関数
function _sphere_detector(z, R, refs, srad, step_size)
    n = size(R,2)
    cand_sets = [Int[] for _ in 1:n] # 候補点のリスト
    x_ml      = zeros(Float64, n)        # ML点ベクトル
    metric    = Inf                      # metric初期値
    cand_list = zeros(Int64, n)
    visited_list = Matrix{Int64}(undef, n, 0)
    metrics = Float64[]
    
    status = 0
    layer  = 1
    while true
        if layer == 0
            if status == 0
                srad *= step_size; layer += 1
            else
                break
            end
        elseif 1<= layer <= n # 探索
            lb, ub = _calc_bounds(z, R, cand_list, refs, srad, layer)
            layer = _search_ml(cand_list, cand_sets, refs, lb, ub, layer)
        elseif layer == n+1 # 
            x_tmp = @views refs[cand_list]
            metric_tmp = _squared_error(z, R, x_tmp)
            push!(metrics, metric_tmp)
            visited_list = hcat(visited_list, cand_list)
            if metric_tmp < metric
                x_ml .= x_tmp
                metric = metric_tmp
            end
            if all(isempty.(cand_sets)) 
                break
            else
                layer = _trace_back(cand_list, cand_sets, refs, layer-1)
            end
            status = 1
        end
    end
    x_ml, metric, metrics, visited_list
end


# 各ステージでの処理
function _search_ml(list, cand_sets, refs, lb, ub, layer)
    n = length(list)
    idx = n - layer + 1
    cand_set = cand_sets[idx]

    # 区間内のシンボルをリストに加える
    ϵ = 1e-3
    for (k,x) in enumerate(refs)
        if (lb - ϵ <= x <= ub + ϵ)
            push!(cand_set, k)
        end
    end

    # 判定
    if isempty(cand_set)
        layer = _trace_back(list, cand_sets, refs, layer-1)
    else
        list[idx] = popfirst!(cand_set)
        layer += 1
    end
    return layer
end


# ２乗距離の計算
function _squared_error(z, R, x)
    e = z - R * x |> r -> dot(r,r)
    return e
end

# 探索した木を逆戻りする
function _trace_back(list, cand_sets, refs, layer)
    n = length(list)

    # 候補点を変更する
    while layer > 0
        idx = n - layer + 1
        cand_set = cand_sets[idx]
        if isempty(cand_set)
            layer -= 1
        else
            list[idx] = popfirst!(cand_set)
            layer += 1
            break
        end
    end
    return layer
end

# 境界条件を計算する
function _calc_bounds(z, R, list, refs, srad, layer)
    n = length(list)
    α = β = 0.0
    i = n - layer + 1 # current index

    # calc alpha
    α_tmp = 0.0
    for k in n:-1:i+1
        temp = z[k]
        for j in k:n
            x_j = refs[list[j]]
            temp -= R[k,j] * x_j
        end
        α_tmp += abs2(temp)
    end
    α = sqrt(srad - α_tmp)

    # calc beta
    β = z[i]
    for j in i+1:n
        x_j = refs[list[j]]
        β -= R[i,j] * x_j
    end
    
    sign_r = sign(R[i,i])
    upper_bound  =  (α + sign_r*β) / abs(R[i,i])
    lower_bound  = (-α + sign_r*β) / abs(R[i,i])

    lower_bound, upper_bound
end


end # module
