get_refs(p::Pam{2}) = pam2_refs
get_refs(p::Pam{4}) = pam4_refs


normfactor(::Pam{2}) = 1.0
normfactor(p::Pam{4}) = ifelse(p.isnorm, sqrt(5.0), 1.0)
normfactor(p::Pam{8}) = ifelse(p.isnorm, sqrt(21.0), 1.0)

normalize!(p::Pam, x) = p.isgrayenc && rdiv!(x, normfactor(p))
denormalize!(p::Pam, x) = p.isgrayenc && rmul!(x, normfactor(p))
