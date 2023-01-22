"""
Qam{M}

example)
    M = 4
    qam = Qam(M)

QAM変復調
    mod(::Qam{M}, inp::AbstractArray{T}) where T
    demod(::Qam{M}, inp::AbstractArray{T}) where T
"""
struct Qam{M}
    bps::Int            # bit/symbol
    isnorm::Bool        # 電力正規化
    isgrayenc::Bool     # Gray符合化
    origin::Symbol      # 起点

    Qam{M}(isnorm, isgrayenc) where M = new{M}(isnorm, isgrayenc, :LL)
    function Qam(M; isnorm=true, isgrayenc=true)
        @assert M % 2 == 0
        new{M}(bitpersymbol(M), isnorm, isgrayenc, :LL)
    end
end

Pam(::Qam{2})  = Pam{2}(false, false)
Pam(q::Qam{4}; isnorm=q.isnorm, isgrayenc=q.isgrayenc)  = Pam{2}(false, isgrayenc)
Pam(q::Qam{16}; isnorm=q.isnorm, isgrayenc=q.isgrayenc) = Pam{4}(false, isgrayenc)

getrefs(q::Qam{2})  = qam2_refs  ./ normfactor(q)
getrefs(q::Qam{4})  = qam4_refs  ./ normfactor(q)
getrefs(q::Qam{16}) = qam16_refs ./ normfactor(q)


"""
シンボル正規化係数の取得
"""
normfactor(q::Qam{2})   = 1.0
normfactor(q::Qam{4})  = ifelse(q.isnorm,  sqrt(2.), 1.0)
normfactor(q::Qam{16}) = ifelse(q.isnorm, sqrt(10.), 1.0)
normfactor(q::Qam{64}) = ifelse(q.isnorm, sqrt(42.), 1.0)

normfactor(::Type{Qam{2}})  = 1.0
normfactor(::Type{Qam{4}},  isnorm)  = ifelse(isnorm,  sqrt(2.), 1.0)
normfactor(::Type{Qam{16}}, isnorm)  = ifelse(isnorm, sqrt(10.), 1.0)
normfactor(::Type{Qam{64}}, isnorm)  = ifelse(isnorm, sqrt(10.), 1.0)

function normfactor(M::Integer, isnorm)
    M == 2  && normfactor(Qam{2}, isnorm)
    M == 4  && normfactor(Qam{4}, isnorm)
    M == 16 && normfactor(Qam{16}, isnorm)
    M == 64 && normfactor(Qam{64}, isnorm)
end
