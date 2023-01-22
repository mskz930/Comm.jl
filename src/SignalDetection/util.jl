# 複素ベクトル(行列)を実ベクトル(実行列)に変換する
to_real(v::AbstractVector{T}) where T<:AbstractFloat = v
to_real(v::AbstractVector{T}) where T<:Complex = vcat(real(v), imag(v))
to_real(V::AbstractMatrix{T}) where T<:AbstractFloat = V
function to_real(V::AbstractMatrix{T}) where T<:Complex
    return [real(V) -imag(V);
            imag(V) real(V)]
end


function soft_output(x_list::AbstractArray, x_metrics::AbstractVectors, S)
    for x in eachrow(x_list)
        
    end
end

function calc_llr_by_list(x, tmp0, tmp1)
    for l = eachindex(idxs)
        if i in S⁺[l]
            tmp1 = ifelse(tmp1>metric[l], metric[l], tmp1)
        else
            tmp0 = ifelse(tmp0>metric[l], metric[l], tmp0)
        end
    end
    return tmp0, tmp1
end