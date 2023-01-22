module Comutils

using Random: randperm, bitrand
import Random: shuffle!
using OffsetArrays: OffsetMatrix
using LinearAlgebra: I
using FFTW: fft!
using RecipesBase



include("eval.jl")

export crandn,
       dftmat, DFT,
       doppler,
       velocity,
       DFTMatrix,
       interleave!,
       deinterleave!,
       decide,
       squared_error,
       show_progress_bar,
       calc_bit_err,
       BER, SER, BLER,
       berplot,
       berplot!,
       BitErrorRates

struct Interleaver
    state::Int64
    n_inputs::Int64
    inds::Vector{Int64}
end

function interleave!(intlv::Interleaver, x)
    x[:] .= x[intlv.inds]
    return x
end

function deinterleave!()
end

function init!(intlv::Interleaver)
    rng = MersenneTwister(intlv.state)
    randperm!(rng, intlv.inds)
    return intlv
end

# 並べ変え順序の初期化
function shuffle!(intlv::Interleaver, x)
    shuffle!(intlv.inds)
    return intlv
end



"""
複素ガウスの正規乱数生成
"""
crandn(dims...) = randn(ComplexF64, dims...)
crandn(::Type{T}, dims...) where T <: Complex = randn(T, dims...)


"""
インターリーブ
"""
function interleave!(input, permind)
    input .= input[permind]
end
interleave(input, permind) = input[permind]

"""
    deinterleave!()

デインタリーバ
"""
function deinterleave!(input, permind)
    input[permind] .= input
    input
end
function deinterleave(input, permind)
    output = similar(input)
    output[permind] = input
    output
end

"""
    doppler(v; fc, cosθ=1.0, scale=:default)

最大ドップラー周波数を計算する

    fc    : center frequency
    v     : velocity
    scale : :default
        default : m/s
        kilo    : convert k/h => m/s
"""
function doppler(v; fc, cosθ=1.0, scale=:default)
    c = 3.0 *  1e8 # 光速
    if scale == :kilo
        v = v * 1e3 / 3600
    end
    fd = (v * fc * cosθ) / c
    return fd
end

# ドップラー周波数から速度を逆算する
function velocity(fd; fc, cosθ=1.0)
    c = 3.0 *  1e8 # 光速
    v = (c * fd) / (fc * cosθ)
    return v
end


"""
get DFT Matrix
"""
function DFT(N, normalize=false)
    F = Array{ComplexF64}(I, N, N)
    F = fft!(F, 1)
    normalize && (F *= 1/sqrt(N))
    F
end

"""
    dftmat(N; normalize=false)
N-point DFT matrix
"""
function dftmat(N, normalize)
    mat = Array{ComplexF64}(I, N, N)
    mat = fft!(mat, 1)
    normalize && (mat *= 1/sqrt(N))
    mat
end
dftmat(N; normalize=false) = dftmat(N, normalize)

"""
DFTMatrix
    W = DFTMatrix(N; normlize=false)
    W[k,n] means exp(im*2pi*k*n/N)
"""
function DFTMatrix(N; normalize=false, T=ComplexF64)
    dftm = dftmat(N, normalize)
    OffsetMatrix(dftm, 0:N-1, 0:N-1)
end


"""
尤度比から判定ビットを取得
"""
decide(x) = x .> 0


"""
誤りビット数計算
"""
calc_bit_err(x, y) = sum(x .!= y)
calc_bit_err(x, y, len) = @views sum(x[1:len] .!= y[1:lem])
function calc_bit_err(x, y, len, stride)
    n_err = 0
    i = j = 1
    while i <= length(y)
        n_err += x[i] == y[j]
        i += 1
        j += stride
    end
    n_err
end

"""
プログレスバーを表示する
"""
function show_progress_bar(n, k=4; progress, digits=3, option="")
    print("process: [" * "-"^(k*n) * " "^(10k-k*n) *"] " *
          "$(round(progress*100, digits=digits))%" * option)
    print("\r")
end


# BERプロット用の構造体とレシピ
struct BER end
struct BLER end
struct SER end
"""
plot(BER(), data)

where the data must be Dict{String,Any} type which has keys{:range, :values, :measure}.
"""
@recipe function f(::BER, data::AbstractDict)
    seriestype := :path
    yscale := :log10
    grid   := true
    markershape --> :auto
    markersize  --> 8
    xguide := "$(data[:measure]) [dB]"
    yguide := "Average BER"
    data[:range], data[:values]
end
@recipe function f(::BLER, data::AbstractDict)
    seriestype := :path
    yscale := :log10
    grid   := true
    markershape --> :auto
    markersize  --> 8
    xguide := "$(data[:measure]) [dB]"
    yguide := "Average BLER"
    data[:range], data[:values]
end



end # module
