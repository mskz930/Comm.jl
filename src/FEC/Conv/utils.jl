
const MINF = 1e20 # 対数の最小値定義

# log_sum_exp
log_sum_exp(x::AbstractVector) = log(sum(exp.(x)))

# jacobian logarithm
jacobian_logarithm(x::Real, y::Real) = max(x,y) + log(1 + exp(-abs(x-y)))
jacobian_logarithm(x::Real, y::Real, z...) = jacobian_logarithm(jacobian_logarithm(x, y), z...)


# Look-Up-Table
struct LookUpTable
    level::Int64
    refs::Vector{Float64}
    vals::Vector{Float64}
    function LookUpTable(;lev=8) # level: 量子化ビット数
        f = x -> log(1 + exp(-x))     # 対象となる関数
        g = y -> -log(exp(y) - 1)     # 逆関数
        step = f(0.0) - f(5.0) / lev
        vals = f.(step .+ step .* [0:lev-1;]) # 量子化(実数値)
        refs = abs.(g.(vals))                 # 参照点
        new(lev, refs, vals)
    end
end

# look-up-tableによる近似log-MAP計算
function (lut::LookUpTable)(x::Real, y::Real)
    z = abs(x-y)
    idx = findlast(lut.refs .> z)
    if isnothing(idx)
        return max(x,y)
    else
        return max(x,y) + lut.vals[idx]
    end
end
