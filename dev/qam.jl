
abstract type Nomralized end
abstract type Unnormalized end

struct Qam{M,T}
  m::Int64
  normalize::Bool
end

normalize(q::Qam{M,Normalize}, x) = normalize(x, nfactor(q))
normalize(q::Qam{M,Unnormalize}, x) = x



