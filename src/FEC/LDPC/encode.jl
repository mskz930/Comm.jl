# encode.jl

# LDPC符号化(function)
# Intのビット数が足りないと結果がオーバーフローしてマイナスになるので注意！
encode(b::AbstractVector, G::AbstractMatrix) = (G * b) .% 2
function encode!(y::T, x::T, G::AbstractMatrix) where {T <: AbstractVector}
  c .= (G * b) .% 2
  return c
end

"""
Encoder: Encoder(G)
    G: Generator Matrix
  LDPC符号化オブジェクト
"""
struct Encoder <: AbstractEncoder
    G::SparseMatrixCSC
end
show(io::IO, e::Encoder) = print(io, "$(typeof(e)): $(size(e)) LDPC Encoder")
show(io::IO, m::MIME"text/plain", enc::Encoder) = print(io, "$(typeof(enc)): $(size(enc)) LDPC Encoder")

function Encoder(G::AbstractMatrix)
  @assert maximum(G) <= 1
  Encoder(sparse(G))
end
# LDPC符号化(functor)
function (enc::Encoder)(inp)
  inp = sparse(inp)
  return Array{Bool}((enc.G * inp) .% 2)
end

code_rate(e::Encoder) = size(e.G,1) // size(e.G,2)
size(e::Encoder) = size(e.G)




