# 発展的なチャネル推定法の関数


# argmax H(i) {Q(H(i)|H(i-1))}
function em_update!(H_new, H_old, x̂, y, β)
    # 平均データ尤度のHに関する極値 ∂/∂H = 0
    # H_new .= H_old
    return
end
# EMアルゴリズムによるCFR係数の更新
function cfr_update_by_em(y, H, β, x̂=nothing, n_iter=1)
    H_old = H
    H_new = similar(H)
    if isnothing(x̂)
        x̂ = Array{eltype(y)}(size(H,2))
        mmse!(x̂, y, H, 1/β)
    else
        mmse!(x̂, y, H, 1/β, x_pri=copy(x̂))
    end
    iter = 1
    while iter < n_iter
        H_new = similar(H_old)
        _update!(H_new, H_old, x̂, y, β)
        iter += 1
    end
    return H_new
end
