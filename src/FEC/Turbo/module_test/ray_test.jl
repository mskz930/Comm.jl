
module Test
include("../Turbo.jl")
using Random: randperm
nbits = 4096
rate = 1//2
turbo = Turbo.TurboCode(4, [13 15 17], 13, bitlen=nbits, iterations=10, rate=rate)
ncbits = Int(turbo.bitlen/rate) + turbo.tailbits
permind = randperm(ncbits)

function main(ebnodBs)
    ebnodBs = Vector(ebnodBs)
    ber_results = zeros(size(ebnodBs))
    for i in eachindex(ebnodBs)
        ber = bercalc(ebnodBs[i])
        @show ebnodBs[i], ber
        if ber > 0
            ber_results[i] = ber
        else
            break
        end
    end
    ber_results
end

function bercalc(ebnodB)
    m_ebnodB = ebnodB + 10*log10(turbo.rate)
    N0 = 10^(-m_ebnodB/10)
    numbits=0; errors=0
    @time while numbits<4e6 && errors<1e5
        tx_bit = rand(nbits) .> 0.5 # 送信ビット
        enc_bit = Turbo.turboenc(turbo, tx_bit)
        enc_bit = enc_bit[permind]
        tx_symbol = 2.0 .* enc_bit .- 1.0
        h = abs.(randn(ComplexF64, length(tx_symbol)))
        # h = abs.(randn(ComplexF64, 1, 1024)) # チャネル係数
        # h = append!(repeat(h,64,1)[:], abs.(randn(ComplexF64, 12)))
        rx_symbol = h.*tx_symbol
        rx_symbol .+= sqrt(N0/2) .* randn(size(tx_symbol))
        Lch = 4.0*(rx_symbol.*h)./N0
        Lch[permind] = Lch
        Lapp = Turbo.turbodec(turbo, Lch)
        rx_bit = Lapp .> 0
        error = sum(tx_bit .!== rx_bit)
        errors += error
        numbits += nbits
    end
    errors/numbits
end

function test(ebnodB)
    N0 = 10^(-ebnodB/10) # N0
    tx_bit = rand(nbits) .> 0.5 # 送信ビット
    enc_bit = Turbo.turboenc(turbo, tx_bit)
    enc_bit = view(enc_bit, permind)
    tx_symbol = 2.0 .* enc_bit .- 1.0
    rx_symbol = tx_symbol .+ sqrt(N0/2) .* randn(size(tx_symbol))
    Lch = 4.0*rx_symbol./N0
    @time Lapp = Turbo.turbodec(turbo, Lch)
end

end
