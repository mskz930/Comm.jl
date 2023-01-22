# awgn_test.jl
using Commun
using Random
Random.seed!(1234) # seed値

G, H = LDPC.load(rows=1024, cols=2048)
const codelen, bitlen = size(G)
const rate = bitlen/codelen
const encoder = LDPC.Encoder(G)
const decoder = LDPC.Decoder(H)
const MAX_BITS = 4e6
const MAX_ERR_BITS = 1e5


function main(;EbN0s)
    qam = Qam(4)
    for EbN0 in EbN0s
        print("EbN0: $(EbN0), ")
        BER = routine(EbN0, qam)
        print("BER: $BER\n")
    end
end

function test(EbN0)
    global rate, bitlen, encoder, decoder
    n_bit = bitlen
    qam = Qam(4)
    SNR = EbN0 + 10*log10(qam.bps * rate)
    N0  = 10^(-SNR/10)
    transceiver(n_bit, N0, qam, encoder, decoder)
end

function routine(EbN0, qam)
    global rate, bitlen, codelen, encoder, decoder, MAX_BITS, MAX_ERR_BITS
    SNR  = EbN0 + 10*log10(qam.bps * rate)
    N0   = 10^(-SNR/10)
    n_bit  = bitlen
    total_bit = total_err_bit = 0
    while total_bit < MAX_BITS && total_err_bit < MAX_ERR_BITS
        n_err_bit = transceiver(n_bit, N0, qam, encoder, decoder)
        total_err_bit += n_err_bit
        total_bit += n_bit
    end
    return total_err_bit/total_bit
end

function transceiver(n_bit, N0, qam, encoder, decoder)
    tx_bits  = randbit(n_bit)
    enc_bits = encoder(tx_bits)
    tx_data  = QAM.mod(qam, enc_bits)
    rx_data = awgn(tx_data, N0)
    rx_LLRs = QAM.demod(qam, rx_data, N0)
    dec_LLRs = decoder(rx_LLRs)
    rx_bits = @views dec_LLRs[1:n_bit] .> 0
    return sum(tx_bits .!= rx_bits)
end
