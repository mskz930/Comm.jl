# common.jl

"""
    randbit(dims...)
ランダムビット列の生成
"""
randbit(dims...) = bitrand(dims...) |> Array

# methods
bitpersymbol(m::Int) = ndigits(m - 1, base = 2)
bps = bitpersymbol


# # PSK型
# struct PskMod
#   M::Int16
#   bps::Int64     # bit/symbol
#   norm::Bool     # 正規化
#   gray::Bool     # Gray符合化
#   ϕ0::Float64    # 初期位相
#   function Psk(M; gray = true, ϕ0 = 0.0)
#     bps = Int(log2(M))
#     Psk(M, bps, false, gray, ϕ0)
#   end
# end


"""
SoftSymbol型
"""
struct SoftSymbol{T<:Number}
  mean::T
  var::Float64
  function SoftSymbol(m, v)
    @assert 0.0 <= v <= 1.0
    new{eltype(m)}(m, v)
  end
  SoftSymbol{T}(m, v) where {T} = new{T}(m, v)
end

Base.show(io::IO, s::SoftSymbol) = print(io, """mean: $(s.mean), :var $(s.var)""")
Base.show(io::IO, m::MIME"text/plain", s::SoftSymbol) = print(io, """$(s.mean)""")

# ソフトシンボル
getmean(x::SoftSymbol{T}) where {T} = getfield(x, :mean)
getvar(x::SoftSymbol{T}) where {T} = getfield(x, :var)

