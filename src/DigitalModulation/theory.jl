using SpecialFunctions: erfc


# QAM変調のビット誤り率理論値
function theory_ber(::Type{Qam{M}}; rng, ratio=:SNR, channel=:AWGN) where M
    ratio==:SNR && (rng = rng .+ 10*log10(log2(M)))
    rng = @. 10^(rng/10)

    m = sqrt(M)
    a = 1/log2(m)*((m-1)/m);
    b = 1/2*3/(m^2-1)

    # BER関数取得: xはEs平均シンボル電力を表す。
    if channel==:AWGN || channel==:awgn
        @. a*erfc(sqrt(b*rng)) # BER関数
    elseif channel==:Rayleigh || channel==:ray
        @. a*(1.0-sqrt(b*rng/(1+b*rng)))
    end
end

# PSK変調のビット誤り率理論値
function theory_ber(::Type{Pam{M}}; rng, ratio=:SNR, channel=:AWGN) where M
    @assert M <= 4
    rng = scale==:SNR ? rng : rng .+ 10*log10(log2(M))  # 電力補正
    rng = 10 .^ (rng/10)                               # dB -> Linear
    m = sqrt(M)
    if M == 2
        a = 1/2; b = 1.0
    else
        a = 1/2; b = 1/2
    end
    # BER関数取得: xはEs平均シンボル電力を表す。
    if channel==:AWGN || channel==:awgn
        @. a*erfc(sqrt(b*x)) # BER関数
    elseif channel==:Rayleigh || channel==:ray
        @. a*(1.0-sqrt(b*x/(1+b*x)))
    end
end

#=
function BER(::Type{QAM}, rng; M, ratio=:SNR, channel=:AWGN)
    rng = ratio==:SNR ? rng : rng .+ 10*log10(log2(M))  # 電力補正
    xlin = 10 .^ (rng/10) # dB -> Linear
    a = 1/log2(M)*((M-1)/M);
    b = 3/(m^2-1)
end
=#


"""
シンボル誤り率の理論値
"""
function theory_ser(::Type{Qam{M}}; rng, channel=:awgn) where M
    rng = collect(rng) # 入力
    x = 10 .^(rng/10)  # Linear scale

    if channel==:awgn || channel==:AWGN

    elseif channel==:ray || channel==:Rayleigh
        @. 2*((sqrt(M)-1)/sqrt(M))*(1-sqrt(3x/(2(M-1)+3x))) - ((sqrt(M)-1)/sqrt(M))^2*(1-sqrt(3x/(2(M-1)+3x))*(4/pi*atan(sqrt((2(M-1)+3x)/3x))))
    end
end
