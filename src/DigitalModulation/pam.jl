
using Binary: bin2dec, grayenc
import Base: mod

struct Natural end
struct Gray end

bitmap(x) = 2.0x - 1.0

"""
Pam{M}
"""
struct Pam{M}
  bps::Int
  normalize::Bool
  grayenc::Bool
  scaling_factor::Float64
end
Pam(M; normalize = true, grayenc = true) = Pam{M}(bitpersymbol(M), normalize, grayenc, ifelse(normalize, normfactor(Pam{M}), 1.0))
function Pam(M::Integer, normalize::Bool = true, grayenc::Bool = true)
  Pam{M}(bitpersymbol(M), normalize, grayenc, ifelse(normalize, normfactor(Pam{M}), 1.0))
end

Base.show(io::IO, pam::Pam{M}) where {M} = print(io,
  "(bps = $(pam.bps), " *
  "normalize = $(pam.normalize), " *
  "grayenc = $(pam.grayenc))"
)


Base.show(io::IO, m::MIME"text/plain", pam::Pam{M}) where {M} = print(
  io,
  """$(typeof(pam)):
  bps       : $(pam.bps)
  normalize : $(pam.normalize)
  grayenc   : $(pam.grayenc)
"""
)


const Pam2_refs = [-1.0, 1.0]
const Pam4_refs = [-3.0, -1.0, 1.0, 3.0]
const Pam8_refs = [-7.0, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0]

get_refs(p::Pam{2}) = Pam2_refs
get_refs(p::Pam{4}) = Pam4_refs

normfactor(p::Pam{2}) = 1.0
normfactor(p::Pam{4}) = ifelse(p.normalize, sqrt(5.0), 1.0)
normfactor(p::Pam{8}) = ifelse(p.normalize, sqrt(21.0), 1.0)
normfactor(::Type{Pam{2}}) = 1.0
normfactor(::Type{Pam{4}}) = sqrt(5.0)
normfactor(::Type{Pam{8}}) = sqrt(21.0)
normalize!(p::Pam, x) = p.grayenc && rdiv!(x, normfactor(p))
denormalize!(p::Pam, x) = p.grayenc && rmul!(x, normfactor(p))


function scaling(p::Pam{M}, m, v) where {M}
  (!p.normalize || M == 2) && return m, v
  m / p.scaling_factor, v / (p.scaling_factor^2)
end
function rescaling(p::Pam{M}, y, σ) where {M}
  (!p.normalize || M == 2) && return y, σ
  y * p.scaling_factor, σ * p.scaling_factor^2
end


"""
Pam変調シンボルの乱数生成
"""
function Base.rand(pam::Pam{M}, size::Integer...) where {M}
  makerefs = getrefs(pam) ./ normfactor(pam)
  rand(refs, size...)
end



# Pamシンボルマップ関数
smap(::Type{Pam{2}}, x::Bool, args...) where {M} = 2x - 1
function smap(::Type{Pam{M}}, b::AbstractVector{T}) where {M,T<:Bool}
  s = 0.0
  sign = 1.0
  val = 0.0
  a = 1 << (length(b) - 1) # 係数
  for l = 1:length(b) # l-th bth
    val = 2 * b[l] - 1
    s += sign * a * val
    sign *= ifelse(grayenc, -val, 1.0)
    a >>= 1
  end
  s
end
function smap(::Type{Pam{M}}, m::Integer, k::Integer) where {M}
  s = 0.0
  a = 1
  for i = 1:k
    s += a * bitmap((m & a) >> (i - 1))
    # @show s
    a <<= 1
  end
  s
end
function smap(::Type{Pam{M}}, λs::Union{T,AbstractArray{T}}, grayenc::Bool = true, normalize::Bool = false) where {M,T<:AbstractFloat}
  a = 1 << (length(λs) - 1) # 係数
  c = 1.0 # 符号
  mu = var = 0.0    # 軟値ビット
  for λ in λs
    s = tanh(λ / 2)
    mu += c * a * s
    var += abs2(a) * (1.0 - abs2(s))
    a >>= 1
    c = ifelse(grayenc, -c * sign(s), c)
    # @show c
  end
  mu, var
end

smap(p::Pam{M}, x::Integer) where {M} = smap(Pam{M}, x, p.bps)
function smap(pam::Pam{M}, λs::Union{Real,AbstractVector{<:Real}}) where {M}
  m, v = smap(typeof(pam), λs, pam.grayenc)
  # return scaling(pam, m, v)
end

# Pam demapping
function sdemap!(::Type{Pam{M}}, x::AbstractVector{T}, y::Real, isgrayenc) where {M,T<:Bool}
  x[1] = y >= 0.0 # MSB
  length(x) == 1 && return
  c = 1 << (length(x) - 1)
  th = ifelse(y > 0, 1.0 * c, -1.0 * c) # しきい値
  for n = 2:length(x)
    x[n] = y >= th
    isgrayenc && (x[n] ⊻= x[n-1])
    c >>= 1
    th += ifelse(y >= th, c, -c)
  end
end
sdemap(pam::Pam{M}, y::Real, σ::Real) where {M} = _approxllr(Pam{M}, rescaling(pam, y, σ)...)

# 参照シンボルを生成
function makerefs(p::Pam{M}) where {M}
  refs = Array{Float64}(undef, M)
  η = normfactor(p)
  for i = 0:M-1
    p.grayenc && (m = grayenc(i))
    refs[i+1] = map(Pam(M), m, p.bps) / η
  end
  return refs
end

"""
    mod(input::AbstractVector{T}; M)

==========
Pam変調
    2-Pam : {0, 1} -> {-1, +1}
    4-Pam : {00,01,10,11} -> {-3, -1, +1, +3}
    8-Pam : {000,001,010,011,100,101,110,111} -> {-7, -5, -3, -1, +1, +3, +5, +7}

"""
function mod(::Pam{2}, inp::AbstractArray{T}) where {T<:Integer}
  return bitmap.(inp)
end
function mod(pam::Pam{M}, inp::AbstractArray{T}) where {M,T<:Bool}
  if size(inp, 1) !== pam.bps # reshape
    inp = reshape(inp, pam.bps, :)
  end
  out = zeros(Float64, size(inp, 2))
  η = normfactor(pam)

  for n in axes(inp, 2)
    m = @views bin2dec(inp[:, n]) |> grayenc
    out[n] = smap(typeof(pam), m, pam.bps) / η
  end
  out
end
function mod(pam::Pam{M}, inp::AbstractArray{T}) where {M, T<:Integer}
  inp = vec(inp)
  out = zeros(Float64, length(inp))
  η = normfactor(pam)

  for n in eachindex(inp)
    m = pam.grayenc ? grayenc(inp[n], pam.bps) : inp[n]
    out[n] = smap(typeof(pam), m, pam.bps) / η
  end
  out
end

function mod(pam::Pam{M}, λs::AbstractArray{T}) where {M,T<:AbstractFloat}
  ms = Vector{Float64}(undef, length(λs))
  vs = Vector{Float64}(undef, length(λs))
  mod!(pam, ms, vs, λs)
end

# Pam復調
function demod!(::Pam{2}, λs::AbstractArray{U}, ys::AbstractArray{T}, σs::AbstractArray{U}) where {T,U}
  for n = 1:length(ys)
    y = ys[n]
    σ = ifelse(length(σs) == 1, σs[1], σs[n])
    λs[n] = _approxllr(Pam{2}, y, σ)
  end
  vec(λs)
end
demod(Pam::Pam{2}, rs) = rs .> 0
function demod(pam::Pam{M}, rs::AbstractArray{T}) where {M,T<:Number}
  bs = Vector{Bool}(undef, pam.bps * length(rs))
  l = pam.bps
  for n in eachindex(rs)
    @views demap!(typeof(pam), bs[(n-1)*l+1:n*l], rs[n], pam.grayenc)
  end
  bs
end
function demod(pam::Pam{M}, rs::AbstractArray{T}, σs::AbstractArray{T}) where {M,T<:Number}
  @assert pam.grayenc
  λs = Vector{Float64}(undef, pam.bps * length(rs))
  l = pam.bps
  for n in eachindex(rs)
    λs[(n-1)*l+1:n*l] .= _approxllr(typeof(pam), rs[n], σs[n])
  end
  λs
end




function slice(::Type{Pam{2}}, r::T) where {T<:Real}
  ifelse(r > 0, 1.0, -1.0)
end
function slice(::Type{Pam{4}}, r::T) where {T<:Real}
  if r >= 0
    if r >= 2
      return 3.0
    else
      return 1.0
    end
  else
    if r >= -2
      return -1.0
    else
      return -3.0
    end
  end
end
function slice(::Type{Pam{8}}, r::T) where {T<:Real}
  if r >= 0
    if r >= 4
      if r >= 6
        return 7.0
      else
        return 5.0
      end
    else
      if r >= 2
        return 3.0
      else
        return 1.0
      end
    end
  else
    if r >= -4
      if r >= -2
        return -1.0
      else
        return -3.0
      end
    else
      if r >= -6
        return -5.0
      else
        return -7.0
      end
    end
  end
end



"""
  _approxllr()

Gray符合化PAMシンボルに対する近似LLRを計算

  max log{p(y|x=+1) / p(y|x=-1)}

arguments:
  N0: complex noise power
  
"""
_approxllr(::Type{Pam{2}}, r::Real, N0::Real) = (4 * r) / N0

function _approxllr(::Type{Pam{4}}, r::Real, N0::Real)
  if r >= 0
    if r >= 2
      (8r - 8) / N0, (-4r + 8) / N0 # (3,-1), (1,3)
    else
      4r / N0, (-4r + 8) / N0 # (1,-1), (1,3)
    end
  else
    if r >= -2
      4r / N0, (4r + 8) / N0 # (1,-1),(-1,-3)
    else
      (8r + 8) / N0, (4r + 8) / N0 # (1,-3), (-1,-3)
    end
  end
end

# 8-Pam
#=
function _approxllr_8Pam(r::Real, N0::Real)
    if r >= 0
        if r>=4
            if r>=6
                return (16r-48, 4r+8, ) ./ N0 # (7,-1), (7,3), (5,7)
            else
                return (16r-48, 4r+8, ) ./ N0 # (5,-1), (5,3), (5,7)
            end
        else
            if r>=2
                return (8r-8, -4r+8, ) ./ N0r # (3,1), (3,5), (1,3)
            else
                return (4r, -4r+8, ) ./ N0 # (1,-1), (1,5), (1,3)
            end
        end
    else
        if r>=-4
            if r>=-2
                return (4r, 4r+8, ) ./ N0 # (-1,1), (-1,-5), (-1,-5)
            else
                return (8r+8,4r+8,) ./ N0 # (-3,1), (-3,-5), (-3,-5)
            end
        else
            if r>=-6
                return (8r+8,4r+8, ) ./ N0  # (-5,1), (-3,-5), (-5.-7)
            else
                return (8r+8,4r+8,) ./ N0 # (-7,1), (-3,-7), (-7,-7)
            end
        end
    end
end
=#


# シンボル確率
function prob(Pam::Type{Pam{M}}, λs) where {M}
  bps = Int(log2(M))
  λs = reshape(λs, bps, :)
  bits = collec(0:M-1)
  Pam.isgrayenc && grayenc!.(bits)
  bitsmat = hcat(dec2bin.(bits, len = bps)...)
  x = llrs' .* map(mapping, bits)
  return exp.(x) ./ sum(exp.(x), dims = 1)
end
