
"""
    zero_padding(v, n)

ゼロパディングする
"""
zero_padding(v::AbstractVector{T}, n) where T = [v; zeros(T, n)]
function zero_padding(A::AbstractArray{T,N}, n; dims=1) where {T,N}
    sizes = collect(size(A))
    sizes[dims] = n
    return cat(A, zeros(sizes...), dims=dims)
end


"""
    squared_error(y,x)

2乗誤差|y - x|^2ノルムを計算
===
    Arguments:
        y, x : Array
    Returns:
        squared error
"""
function squared_error(y, x)
    @assert size(x)[2:end]==size(y)[2:end]
    y = reshape(y, size(y,1), :)
    x = reshape(x, size(x,1), :)

    r = size(x,1) - size(y,1)
    if r > 0
        y = zero_padding(y, r, dims=1)
    elseif r < 0
        x = zero_padding(x, -r, dims=1)
    end
    sqerr = 0.0
    for n in axes(x,2)
        sqerr += @views squared_error(x[:,n], y[:,n])
    end
    return sqerr
end

function squared_error(v1::AbstractVector, v2::AbstractVector)
    l1, l2 = length(v1), length(v2)
    if l1 != l2
        r = l1 - l2
        if r > 0
            v2 = zero_padding(v2, r)
        else
            v1 = zero_padding(v1, -r)
        end
    end

    sqerr = 0.0
    for i in eachindex(v1)
        sqerr += abs2(v1[i]-v2[i])
    end
    sqerr
end
