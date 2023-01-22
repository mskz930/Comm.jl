# sphere_decoding.jl
using LinearAlgebra: qr, dot


"""
    sphere_detector(y, H, modparams, radius_init, dr=2) where T

MLdetection by sphere decoding algorithm
    - PAM/QAM modulated symbol is assummed.
"""
function sphere_detector(y::AbstractVector{T}, H::AbstractMatrix{T}, refs, slicer; step_size=2.0) where T <: Complex
    y = [real(y); imag(y)]
    H = [real(H) -imag(H);
            imag(H) real(H)]
    Q, R = qr(H)
    sphere_detector(y, Q, R, refs, sliceer, step_size)
end
function sphere_detector(y::AbstractVector{T}, H::AbstractMatrix{T}, refs, slicer; step_size=2.0) where T <: Real
    Q, R = qr(H)
    sphere_detector(y, Q, R, refs, slicer, step_size)
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
    cand_lists = [Int[] for _ in 1:n] # 候補点のリスト
    x_tmp  = zeros(Float64, n)        # 候補点ベクトル
    x_ml   = zeros(Float64, n)        # ML点ベクトル
    metric = Inf                      # metric初期値
    
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
            lb, ub = _calc_bounds(z, R, x_tmp, srad, layer)
            layer = _search_ml(x_tmp, cand_lists, refs, lb, ub, layer)
        elseif layer == n+1 # 
            metric_tmp = _squared_error(z, R, x_tmp)
            if metric_tmp < metric
                x_ml .= x_tmp
                metric = metric_tmp
            end
            if all(isempty.(cand_lists)) 
                break
            else
                layer = _back_tracking(x_tmp, cand_lists, refs, layer-1)
            end
            status = 1
        end
    end
    x_ml, metric
end