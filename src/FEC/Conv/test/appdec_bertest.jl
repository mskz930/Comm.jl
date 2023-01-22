using Random, Comm
using UnPack
using .FEC.Conv, .Modulator
using .Comutils: interleave!, deinterleave!, calc_bit_err

# BER計算
function set_params(; kwargs...)
    nbits = 1000
    qam = Qam(2)
    trellis = poly2trellis([3 7]); # トレリス
    encoder = inp -> convenc(inp, trellis=trellis)
    decoder = APPDecoder(trellis, numbits=nbits); # MAP復号器パラメータ
    ncbits = numofoutputs(trellis, nbits)

    params = (
        nbits = nbits,
        ncbits = ncbits,
        rate = getrate(trellis),
        mod_t = qam,
        encoder = encoder,
        decoder = decoder
    )
    return params
end

# BERの計算
function main(;EbN0, kwargs...)
    BER = Float64[]
    params = set_params(; kwargs...)
    for ebno in EbN0
        println("EbN0[dB]: $ebno")
        ber = BERcalc(ebno, params)
        if ber > 0
            println("BER: $BER")
            push!(BER, ber)
        else
            break
        end
    end
    (
    range = EbN0[eachindex(BER)],
    value = BERs,
    measure = :EbN0
    )
end

function BERcalc(ebno_dB, params)
    @unpack rate, mod_t, nbits, ncbits = params
    snr_dB = ebno_dB + 10log10(rate * mod_t.bps)
    nvar = 10^(-snr_dB/10)

    total_bits = total_err_bits = iter = 0
    while (total_bits < 4e6) && total_err_bits < 1e4
        n_err_bits = transceiver(; nvar, params...)
        total_bits += nbits
        total_err_bits += n_err_bits
        iter += 1
        if iter % 1000 == 0
            @show iter, total_err_bits, total_bits
        end
    end
    total_err_bits / total_bits
end

# 動作テスト
function test(;SNRdB=10, kwargs...)
    nvar=10^(-SNRdB/10)
    params = set_params(kwargs...)
    transceiver(;nvar, params...)
end

function transceiver(; nvar, nbits, ncbits, rate, encoder, mod_t, decoder)
    perminds = randperm(ncbits)

    tx_bits = randbit(nbits) # バイナリデータ生成
    enc_bits = encoder(tx_bits) # 符号化ビット列
    enc_bits = interleave!(enc_bits, perminds)
    tx_data = mod(mod_t, enc_bits) # バイナリ => シンボル

    rx_data = awgn(tx_data, nvar); # 受信シンボル
    demod_LLRs = demod(mod_t, rx_data, nvar) # 受信LLR
    deinterleave!(demod_LLRs, perminds)
    dec_bits = decoder(demod_LLRs) # 最大事後確率復号
    rx_bits = dec_bits[1:nbits] .> 0
    return calc_bit_err(tx_bits, rx_bits)
end
