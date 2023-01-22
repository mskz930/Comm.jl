# models.jl
module PDP

using OffsetArrays
using StatsBase: sample

export nzinds, rms


"""
    equiuniform(;interval, L)
一様パスモデル
"""
function equiuniform(;L, space=1)
    inds = 1:space:L+1
    coefs = zeros(Float64,last(inds))
    coefs[inds] .= 1.0
    coefs ./= sum(coefs)
    return coefs
end


"""
    random_uniform(k, )

    k : number of non-zero path
    L : maximum channel length

"""
function random_uniform(k, L)
    @assert k > 0
    coefs = zeros(L+1)
    coefs[1] = 1.0
    coefs[sample(2:L+1, k-1, replace=false)] .= 1.0
    coefs ./= sum(coefs) # 正規化
    coefs
end

"""

指数減衰チャネルモデルの生成
    arguments:
        L : 最大遅延パス数(チャネル次数)
"""
function exponent(;L, ratio, space=1)
    L==0 && return [1.0]
    
    # パラメータ
    A = 10^(-ratio/10)   # dB => linear
    inds = 1:space:L+1   # パス位置
    L = last(inds)-1     # チャネル次数
    coefs = zeros(L+1)

    # プロファイルの作成
    τ_m = length(inds)-1                          # 遅延パス数
    σ_τ = τ_d = -τ_m/log(A)                         # 最大遅延スプレッド
    p0 = (1-exp(-(τ_m+1)/σ_τ))/(1-exp(-1/σ_τ))      # 正規化定数
    coefs[inds] .= 1/p0 * exp.(-(0:τ_m)/σ_τ)       # パス係数

    return coefs[1:maximum(inds)]
end

# 平均遅延時間の計算
function mean(pdp)
    τmax = length(pdp)-1
    μ = 0.0
    for l in 0:τmax
        μ += l * pdp[l+1]
    end
    μ/τmax
end

# 標準偏差
function rms(pdp)
    μ = mean(pdp) # 平均
    σ2 = 0.0 # 分散
    τmax = length(pdp)-1
    for l in 0:τmax
        σ2 += (l-μ)^2 * pdp[l+1]
    end
    return sqrt(σ2)
end

# 非ゼロなパスのindexを返す
nzinds(coefs) = findall(x -> !iszero(x), coefs)




end # module
