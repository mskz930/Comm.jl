import Base: map, mod

export bpskmod, qpskmod, bpskdemod, qpskdemod

"""
  bpskmod()

BPSKmodulation

arguments: 
  x : {0,1}
returns:
  BSPK modulated symbol

"""
bpskmod(x::Integer) = 2.0 * (x % 2) - 1.0
bpskmod(x::AbstractArray) = bpskmod.(x)

"""
  qpskmod(x::Integer)

arguments:
  x : {0,1,2,3}

returns:
  QPSK modulated symbol
"""
function qpskmod(x::Integer)
  s1 = 2.0 * (x & 1) - 1.0
  s2 = 2.0 * (x & 2 >> 1) - 1.0
  return s1 + im * s2
end
qpskmod(x::AbstractArray) = qpskmod.(x)

bpskdemod(ys::AbstractArray{T}) where T <: Number = bpskdemod.(ys)

function qpskdemod(ys::AbstractArray)
  bs = Matrix{Bool}(undef, 2, length(ys))
  for i = 1:size(bs, 2)
    bs[:, i] = qpskdemod(ys[i])
  end
  bs
end

struct PSK{M} end

# BPSK = PSK{2}
# QPSK = PSK{4}

