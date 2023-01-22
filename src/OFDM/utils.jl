# OFDM/src/common.jl

# 配列をサイズnのサブ配列に分割
function _split(v::AbstractVector{T}, n::Integer) where T
    vs = [v[i:n:end] for i in 1:n]
    return vs
end

# elmsの配列の値の要素数のサブベクトルに分割する
# [elems[1], elems[2], elems[3]...]
function _split(x::AbstractArray{T}, elms::Vector{<:Integer}) where T
    length(elms) == 1 && return [x]
    length(x) != sum(elms) && error("elmsの総数は配列の長さに一致する必要があります")
    arr = Vector{T}[]; sizehint!(arr, length(elms))
    n = 1
    for i in eachindex(elms)
        subarr = T[]
        for j in 1:elms[i]
            push!(subarr, x[n])
            n += 1
        end
        push!(arr, subarr)
    end
    return arr
end


# ベクトルを複数のサブベクトルに並行に分割する
# [x1,x2,x3,...] => [[x1,...], [x2,...], [x3...]]
function _to_substreams(v::Vector{T}, n::Integer) where T
    [v[i:n:length(v)] for i in 1:n]
end


# Serial to Parallel
function serial_to_parallel(v::AbstractVector{T}, n, pad=:zero) where T
    len = length(v)
    len == 0 && return v
    rem = length(v) % n # 割り切れるか?
    rem > 0 && _padding!(v, n-rem, pad)
    return reshape(v, n, div(length(v),n))
end

# データパディング
function _padding!(v::AbstractArray{T}, k, by) where T
    if by ==:zero
        append!(v, zeros(T,k))
    elseif by == :elm
        append!(v, rand(v,k))
    end
    v
end

"""
Parallel to Serial
"""
parallel_to_serial(v::AbstractArray{T,2}, divs) where T = vec(v)
parallel_to_serial(v::AbstractArray{T,3}, divs) where T = reshape(v, :, size(v,3))



"""
時間領域受信信号列をstackしてシンボルごとに並べた行列を作成する
"""
function stack(rx_sig, n_fft, n_gi, nrx)
    rx_sig = transpose(rx_sig)
    rx_sig = reshape(rx_sig, (n_fft+n_gi)*nrx, :) # (Nrx*n_fft, Nsym)
    return rx_sig
end

"""
チャネル行列の生成

    Arguments:
        h : channel impulse response vector
        M : symbol block size
        N : FFT length
        G : GI length
"""
function make_channel_matrix(h, M, N, G)
    L, nrx, ntx = size(h,1), size(h,2), size(h,3)
    R = L-1-G
    hc = permutedims(h, [2,1,3])
    hc = reshape(hc, nrx*L, :)
    hc = zero_padding(hc, M-nrx*L, dims=1) # column vector
    hr = permutedims(h, [2,3,1])
    hr = reshape(hr, nrx, :)
    hr = zero_padding(hr, N-ntx*L, dims=2)
    hr = circshift(hr, (0, -nrx))
    hr = My.Linalg.reverse(hr, shape=(nrx, ntx), dir=:row)
    H0 = My.Linalg.toeplitz(hc, hr; shape=(nrx, ntx))

    H10 = H11 = nothing
    if R > 0
        m, n = nrx, ntx
        H10 = zero(H0); H11 = zero(H0)
        for i = 1:m:(L-1-G)*m
            for j = N-n*(L-1)+1+i-1:n:N-n*G
                H10[i:i+m-1,j+G*n:j+n-1+G*n] .= H0[i:i+m-1,j:j+n-1]
                if G > 0
                    H11[i:i+m-1,j:j+n-1] .= H0[i:i+m-1,j:j+n-1]
                end
                H0[i:i+m-1,j:j+n-1] .= zero(eltype(H0))
            end
        end
    end
    return H0, (H10, H11)
end


function make_toeplitz_matrix(elms::AbstractVector, nrows, ncols)
    A = zeros(eltype(elms), nrows, ncols)
    L = length(elms)

    for j = 1:ncols
        m = ifelse(j+L-1<=nrows, j+L-1, nrows)
        for i = j:m
            A[i,j] = elms[i-j+1]
        end
    end
    return A
end

function make_toeplitz_matrix(elms1::AbstractVector, elms2::AbstractVector)
    A = zeros(eltype(elms), nrows, ncols)
    L = length(elms)

    for j = 1:ncols
        m = ifelse(j+L-1<=nrows, j+L-1, nrows)
        for i = j:m
            A[i,j] = elms[i-j+1]
        end
    end
    return A
end
