"""
    transversal_filter

トランスバーサルフィルタ
"""
function transversal_filter(input::AbstractArray, coeffs::AbstractArray)
    promoted_type = promote_rule(eltype(input),eltype(coeffs)) # 優先型
    output = zeros(promoted_type, size(input,1)+size(coeffs,1)-1) # 出力配列
    transversal_filter!(output, input, coeffs) # 畳み込み
    return output
end

function transversal_filter!(out::AbstractVector{T}, inp::AbstractVector, w::AbstractVector) where T
    N = size(inp,1) # 入力長
    L = size(w,1)   # フィルタ長
    M = size(out,1) # 出力長
    @assert M >= N

    # 畳み込み
    @inbounds for m in 1:N
        if m <= L
            mem = @view inp[m:-1:1]
        elseif L < m <= N
            mem = @view inp[m:-1:m-L+1]
        else
            mem = @view inp[end:-1:m]
        end
        val = zero(T)
        for (mi,wi) in zip(mem, w)
            val += mi * wi
        end
        out[m] += val
    end
    return
end

function transversal_filter!(out::AbstractVector{T}, inp::AbstractVector, w::AbstractMatrix) where T
    N = size(inp,1) # 入力長
    L = size(w,1)   # フィルタ長
    M = size(out,1) # 出力長
    @assert M >= N

    # 畳み込み
    @inbounds for m in 1:N
        if m <= L
            mem = @view inp[m:-1:1]
        elseif L < m <= N
            mem = @view inp[m:-1:m-L+1]
        else
            mem = @view inp[end:-1:m]
        end
        val = zero(T)
        for (mi,wi) in zip(mem, view(w,:,m))
            val += mi * wi
        end
        out[m] += val
    end
    out
end
