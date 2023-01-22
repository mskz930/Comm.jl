

include("../turbo.jl")
function main()
    nbits = 1000
    trellis = turbo.poly2trellis(5, [33 37], 33)
    Turbo = turbo.Turbo(trellis, bitlen=nbits, iterations=6)

    ebnodBs = collect(-5:0.5:0)
    ber_results = zeros(Float64, size(ebnodBs))
    for (i,ebnodB) in enumerate(ebnodBs)
        N0 = 10^(-ebnodB/10)
        errors=0
        numbits=0
        while (numbits < 4e6) && (errors<1e3)
            tx_data = rand(nbits) .> 0.5
            tx_code = turbo.turboenc(Turbo, tx_data)
            tx_symbol = 2.0 .* tx_code .- 1.0
            rx_symbol = tx_symbol .+ N0/2 .* randn(size(tx_code))
            Lch = 4.0*rx_symbol./N0
            Lapp, Lext = turbo.turbodec(Turbo, Lch)
            rx_data = Lapp .> 0
            errors += sum(rx_data .!== tx_data)
            numbits += nbits
            if (numbits)%Int(1e6) == 0
                println(numbits)
            end
        end
        ber = errors/numbits
        @show ebnodB, ber
        ber_results[i] = ber
        if ber <= 0
            break
        end
    end
    ebnodBs, ber_results
end



#=
function main()
    ebnodBs= 0:0.5:5
    ber_results = length(ebnodBs)
    ber_results = bit_error(ebnodBs, ber_results)
end

function bit_error()
end
=#
