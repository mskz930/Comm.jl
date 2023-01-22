# utils.jl

# ベクトルを指定の(step,len)によってベクトル配列に分割する
function _split_into_subgroups(arr::Vector{T}; step, len) where T
    st = [1:step:len;]
    en = [step:step:len;]
    length(st)>length(en) && push!(en, len)
    [filter(x-> i<=x<=j, arr) for (i,j) in zip(st,en)]
end


"""
    zero_padding(arr::AbstractArray{T}; pad_size, dims=1) where T

配列にゼロパディングする
"""
function zero_padding(A::AbstractArray{T}, pad_size; dims=1) where T
    sizes = [size(A)...]
    sizes[dims] = pad_size
    pad = zeros(T, sizes...)
    cat(A, pad, dims=dims)
end

"""

パイロット推定誤差分散
"""
function pilot_estimation_error(pilot::Pilot)
end

# ベクトルの2乗誤差
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
    return sqerr / size(x,2)
end

function squared_error(v1::AbstractVector, v2::AbstractVector)
    @assert length(v1) == length(v2)
    error = 0.0
    for i in eachindex(v1)
        error += abs2(v1[i]-v2[i])
    end
    error/length(v1)
end
