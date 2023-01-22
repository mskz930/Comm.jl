# LDPCモジュールの動作テスト

using Commun.FEC.LDPC
using Commun.Digimod
using Commun.Comch: awgn_channel!
using Random: bitrand


G, H = ldpc_code_loader(rows=256, cols=512, method="PEG")
ldpcenc = LDPCEncoder(G)
ldpcdec = LDPCDecoder(H, dectype=:log_sum_product)

n_bit = 256
modtype = Qam(2)

SNRdB = 10
N0  = 10^(-SNRdB/10)

info_bits  = bitrand(256)
enc_bits   = ldpcenc(info_bits)
tx_data    = mod(modtype, enc_bits) |> real
rx_bits    = awgn_channel!(tx_data, N0)
demod_LLRs = demod(modtype, rx_bits, N0)

Juno.@enter ldpcdec(demod_LLRs)
