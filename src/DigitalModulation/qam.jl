import Base: mod, rand


"""
Qam{M}
"""
struct Qam{M}
  bps::Int           # bit/symbol
  normalize::Bool  # normalization flag
  grayenc::Bool    # gray encoding flag
  scaling_factor::Float64 # scaling factor
end
Qam(M::Integer; normalize = true, grayenc = true) = Qam{M}(bitpersymbol(M), normalize, grayenc, ifelse(normalize, normfactor(Qam{M}), 1.0))

Base.show(io::IO, m::MIME"text/plain", qam::Qam{M}) where M = println(io, """
$(typeof(qam)):
  bps       : $(qam.bps)
  normalize : $(qam.normalize)
  grayenc   : $(qam.grayenc)  
""")

Base.eltype(::Qam{M}) where {M} = ComplexF64


# Qam型からPam型を作成する関数

Pam(::Qam{2}) = Pam{2}(false, false)
Pam(q::Qam{4}) = Pam{2}(false, grayenc)
Pam(q::Qam{16}) = Pam{4}(false, grayenc)




const qam2_ref = [-1.0 + 0.0im, 1.0 + 0.0im]
const qam4_ref = [-1.0 - 1.0im, -1.0 + 1.0im, 1.0 - 1.0im, 1.0 + 1.0im]
const qam16_ref = [
  -3.0 - 3.0im, -3.0 - 1.0im, -3.0 + 1.0im, -3.0 + 3.0im,
  -1.0 - 3.0im, -1.0 - 1.0im, -1.0 + 1.0im, -1.0 + 3.0im,
  1.0 - 3.0im, 1.0 - 1.0im, 1.0 + 1.0im, 1.0 + 3.0im,
  3.0 - 3.0im, 3.0 - 1.0im, 3.0 + 1.0im, 3.0 + 3.0im
]
const qam64_ref = [
  -7.0 - 7.0im, -7.0 - 5.0im, -7.0 - 3.0im, -7.0 - 1.0im, -7.0 + 1.0im, -7.0 + 3.0im, -7.0 + 5.0im, -7.0 + 7.0im,
  -5.0 - 7.0im, -5.0 - 5.0im, -5.0 - 3.0im, -5.0 - 1.0im, -5.0 + 1.0im, -5.0 + 3.0im, -5.0 + 5.0im, -5.0 + 7.0im,
  -3.0 - 7.0im, -3.0 - 5.0im, -3.0 - 3.0im, -3.0 - 1.0im, -3.0 + 1.0im, -3.0 + 3.0im, -3.0 + 5.0im, -3.0 + 7.0im,
  -1.0 - 7.0im, -1.0 - 5.0im, -1.0 - 3.0im, -1.0 - 1.0im, -1.0 + 1.0im, -1.0 + 3.0im, -1.0 + 5.0im, -1.0 + 7.0im,
  1.0 - 7.0im, 1.0 - 5.0im, 1.0 - 3.0im, 1.0 - 1.0im, 1.0 + 1.0im, 1.0 + 3.0im, 1.0 + 5.0im, 1.0 + 7.0im,
  3.0 - 7.0im, 3.0 - 5.0im, 3.0 - 3.0im, 3.0 - 1.0im, 3.0 + 1.0im, 3.0 + 3.0im, 3.0 + 5.0im, 3.0 + 7.0im,
  5.0 - 7.0im, 5.0 - 5.0im, 5.0 - 3.0im, 5.0 - 1.0im, 5.0 + 1.0im, 5.0 + 3.0im, 5.0 + 5.0im, 5.0 + 7.0im,
  7.0 - 7.0im, 7.0 - 5.0im, 7.0 - 3.0im, 7.0 - 1.0im, 7.0 + 1.0im, 7.0 + 3.0im, 7.0 + 5.0im, 7.0 + 7.0im,
]


get_refs(q::Qam{2}) = qam2_refs ./ normfactor(q)
get_refs(q::Qam{4}) = qam4_refs ./ normfactor(q)
get_refs(q::Qam{16}) = qam16_refs ./ normfactor(q)



"""
シンボル正規化係数の取得
"""


normfactor(::Type{Qam{2}}) = 1.0
normfactor(::Type{Qam{4}}) = sqrt(2.0)
normfactor(::Type{Qam{16}}) = sqrt(10.0)
normfactor(::Type{Qam{64}}) = sqrt(42.0)


"""
QAMシンボルの正規化/逆正規化
"""
normalize(qam, x) = qam.grayenc ? (x / normfactor(qam)) : x
denormalize(qam, x) = qam.grayenc ? (x * normfactor(qam)) : x
normalize!(qam, x) = qam.grayenc && rdiv!(x, normfactor(qam))
denormalize!(qam, x) = qam.grayenc && rmul!(x, normfactor(qam))

scaling(qam::Qam{M}, x::Number) where M = x / qam.scaling_factor
scaling(qam::Qam{M}, m::Number, v::Number) where M= (m, v) ./ qam.scaling_factor
rescaling(qam::Qam{M}, x::Number) where M = x * qam.scaling_factor
rescaling(qam::Qam{M}, y::Number, σ::Number) where M = (y*qam.scaling_factor, σ*qam.scaling_factor^2)

#======= QAM　変調 ========#
# ビット列を対応するQam信号点へマッピングする
#=
function bitmap(::Type{QAM{M}}, inp) where {M}
  k = bps(M)
  η = qam.norm ? normfactor(qam) : 1.0
  points = Array{ComplexF64}(undef, length(inp))
  if M == 2
    for i = 0:M-1
      refs[i+1] = pam_mapping(i, 1, false)
    end
  else
    sqrtm = 2^(bps >> 1)
    cnt = 1
    hbps = bps >> 1
    for m1 = 0:sqrtm-1
      for m2 = 0:sqrtm-1
        si = pam_mapping(m1, hbps, qam.grayenc)
        sq = pam_mapping(m2, hbps, qam.grayenc)
        points[cnt] = complex(si, sq) / η
        cnt += 1
      end
    end
  end
  points
end
=#


# M-PAMシンボル列からM-QAMシンボル列を構成する
# 実部虚部を交互に順に並べたもの引数として扱う
function Base.merge(::Type{Qam}, real_input::AbstractArray{T}) where {T<:Real}
  len = length(real_input)
  output = Array{Complex{T}}(undef, div(len, 2))
  for i = 1:2:length(real_input)
    output[div(i, 2)+1] = complex(real_input[i], real_input[i+1])
  end
  output
end


"""
QAM変調
"""
function mod(::Qam{2}, inp::AbstractArray{T}) where {M,T<:Integer}
  out = 2.0 * inp .- 1
  return out
end
function mod(qam::Qam{M}, inp::AbstractArray{T}) where {M,T<:Bool}
  if size(inp, 1) != qam.bps
    inp = reshape(inp, Int(qam.bps), :)
  end
  out = zeros(ComplexF64, size(inp, 2))
  k = qam.bps >> 1
  η = normfactor(qam)
  pam = Pam{M >> k}
  for n in eachindex(out)
    # binary -> integer
    m_i, m_q = @views bin2dec(inp[1:k, n]), bin2dec(inp[k+1:end, n])

    x_i = smap(pam, m_i, k)
    x_q = smap(pam, m_q, k)
    out[n] = complex(x_i, x_q) / η
  end
  out
end

# QAM変調(軟値シンボル生成)
function mod!(qam::Qam{M}, ms::AbstractArray{T}, vs::AbstractArray{U}, λs::AbstractArray{U}) where {M,T,U}
  if size(λs, 1) != qam.bps
    λs = reshape(λs, qam.bps, :)
  end

  η = normfactor(qam)
  l = qam.bps >> 1
  pam = Pam{M >> l}

  for n in axes(λs, 2)
    si, vi = @views map(pam, λs[1:l, n], qam.grayenc)
    sq, vq = @views map(pam, λs[l+1:end, n], qam.grayenc)

    ms[n] = (si + sq * im) / η
    vs[n] = (vi + vq) / η^2
  end
  ms, vs
end
function mod(qam::Qam{M}, λs::AbstractArray{T}) where {M,T<:AbstractFloat}
  λs = reshape(λs, qam.bps, :)
  ms = Array{ComplexF64}(undef, size(λs, 2))
  vs = Array{Float64}(undef, size(λs, 2))
  mod!(qam::Qam{M}, ms, vs, λs)
end

# Qamシンボル乱数の生成
function rand(qam::Qam{M}, args::Integer...) where {M}
  refs = get_refs(qam)
  return rand(refs, args...) # random indices: 1~M
end



# QAMシンボルマップ
function smap(qam::Qam{M}, x::Integer) where {M}
  qam.grayenc && (x = grayenc(x, qam.bps))
  l = (qam.bps >> 1) 
  i = (l << 1) - 1
  q = (M-1) - i
  x_re = smap(Pam{M}, (x&q) >> l, l)
  x_im = smap(Pam{M}, (x&i), l)

  return complex(x_re, x_im)
end

function sdemap(qam::Qam{M}, y::Number, σ::Real) where M
  if qam.normalize
    y, σ = rescaling(qam, y, σ)
  end
  l = qam.bps >> 1
  _approxllr(Pam{M >> l}, real(y), σ), _approxllr(Pam{M >> l}, imag(y), σ)
end

# マップリスト(辞書)を作成する
function make_mapping_table(qam::Qam{M}) where {M}
  maplist = OrderedDict{Int,ComplexF64}()
  bps = bitpersymbol(qam)
  for i = 0:M-1
    maplist[i] = qammap(i, bps, qam.grayenc) / normfactor(qam)
  end
  return maplist
end

# 参照点の生成(バイナリの整数表現とシンボルが一致)
function make_refs(qam::Qam{M}) where {M}
  refs = ComplexF64[]
  for i = 0:M-1
    x = qammap(i, bitpersymbol(qam), qam.grayenc) / normfactor(qam)
    push!(refs, x)
  end
  return refs
end




#======= 
QAM復調 
========#

"""
    demod(qam::Qam{M}, y::AbstractArray{T}) where {T<:Number}

QAMシンボル復調
"""
demod(qam::Qam{2}, ys::AbstractArray{T}) where {T<:Number} = demod(Pam(qam), inp)

function demod(qam::Qam{M}, ys::AbstractArray{T,N}) where {M,T<:Complex,N}
  η = normfactor(qam)
  k = qam.bps >> 1
  pam = Pam{M >> k}
  out = Array{Bool}(undef, bps, length(ys))
  for n in eachindex(inp)
    y = y[n] * η
    # Real/Imag => Binary
    @views demap!(pam, out[1:k, n], real(y), qam.grayenc)     # 実部
    @views demap!(pam, out[k+1:end, n], imag(y), qam.grayenc) # 虚部
  end
  vec(out)
end
function demod(qam::Qam{M}, inp::AbstractArray{T,N}) where {M,N,T<:Number}
  inp .*= normfactor(qam)
  bps = bitpersymbol(qam)
  hbps = bps >> 1
  sqrtm = M >> div(bps, 2)
  pam = Pam(qam)
  out = Array{Bool}(undef, bps, length(inp))
  for n = 1:2:length(inp)
    @views demap!(out[1:hbps, n], inp[n], qam.grayenc) # 実部
    @views demap!(out[hbps+1:end, n], inp[n+1], qam.grayenc) # 虚部
  end
  vec(out)
end

# Soft QAM Demodulation
function demod!(qam::Qam{2}, λs::AbstractArray{T}, ys::AbstractArray{U}, σs::AbstractArray{T}) where {M,T,U}
  demod!(Pam(2), λs, ys, σs)
end
demod(qam::Qam{2}, ys::AbstractArray{T}, N0) where {T<:AbstractFloat} = demod(Pam(2; normalize = false), ys, N0)

function demod!(qam::Qam{M}, λs::AbstractArray{T}, ys::AbstractArray{U}, σs::AbstractArray{T}) where {M,T,U}
  @assert qam.grayenc
  if size(λs, 1) != qam.bps
    λs = reshape(λs, qam.bps, :)
  end

  η = normfactor(qam) # 正規化定数
  l = qam.bps >> 1    # 半分のbps
  pam_t = Pam{M >> l} # subtype
  for n = 1:length(ys)
    r = ys[n]
    σ = ifelse(length(σs) == 1, σs[1], σs[n])
    if qam.normalize
      r *= η
      σ *= η^2
    end
    λs[1:l, n] .= _approxllr(pam_t, real(r), σ)
    λs[l+1:end, n] .= _approxllr(pam_t, imag(r), σ)
  end
  @show λs
  vec(λs)
end
function demod(qam::Qam{M}, ys::AbstractArray, σs) where {M}
  λs = Matrix{Float64}(undef, qam.bps, length(ys)) # 出力配列
  demod!(qam, λs, ys, σs)
end





# QAMソフトシンボルを生成
function SoftSymbol(qam::Qam{M}, LLRvec::AbstractVector{T}) where {M,T<:AbstractFloat}
  bps = bitpersymbol(qam)
  η = normfactor(qam)
  mi, vi = @views soft_pam_mapping(LLRvec[1:div(bps, 2)], qam.grayenc)
  mq, vq = @views soft_pam_mapping(LLRvec[div(bps, 2)+1:bps], qam.grayenc)
  mean = (mi + mq * im) / η
  var = (vi + vq) / abs2(η)
  SoftSymbol{ComplexF64}(mean, var)
end

#============================= 
 Slicer
=============================#

slice(::Type{Qam{2}}, r::Complex) = slice(Pam{2}, real(r))
slice(::Type{Qam{M}}, r::Complex) where {M} = slice(Pam{M}, real(r)) + im * slice(Pam{M}, imag(r))


# ソフトシンボル配列の生成
function _merge_soft_symbols(means::AbstractVector{T}, vars::AbstractVector{U}) where {T,U}
  outputs = Vector{SoftSymbol{T}}(undef, length(means))
  for i in eachindex(means)
    outputs[i] = SoftSymbol{T}(means[i], vars[i])
  end
  return outputs
end

# 各ビットが1に対応するシンボル(10進数表現)の集合を返す
function make_binary_sets(qam::Qam{M}) where {M}
  bps = bitpersymbol(qam)
  bsets = Dict{Int,Set{Int}}()
  decs = collect(0:M-1) |> grayenc!
  decs = grayenc!(decs)
  for i = 1:bps
    digit = 1 << (bps - i)
    set = Set{Int}()
    for d in decs
      (d & digit) > 0 && push!(set, d)
    end
    bsets[i] = set
  end
  return bsets
end
