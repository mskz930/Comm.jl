module Utils


# テプリッツ行列の生成
function toeplitz(h::AbstractVector{T}, m, n) where T
    mat = zeros(T, m, n)
    for i in eachindex(h)
        i > m && break
        for j in 1:size(mat,2)
            j+i-1>m && break
            mat[j+i-1,j] = h[i]
        end
    end
    mat
end





end # module
