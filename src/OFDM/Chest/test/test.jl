using Revise
using Random
using BenchmarkTools
using Commun.OFDM
using Commun.Digimod
using Commun.Comch
using Commun.OFDM.Chest: pilot_symbol_average, cir_estimate

function test(SNRdB)
    global BIT_LENGTH, N_DIMS, ofdm, modtype, ofdmframe

    n_dims = N_DIMS
    bitlen = BIT_LENGTH

    tx_frame = copy(ofdmframe)
    N0 = 10^(-SNRdB/10)

    tx_bits = bitrand(bitlen)
    tx_data = mod(modtype, tx_bits)
    tx_sig  = ofdm_mod(ofdm, tx_data, tx_frame)

    rx_sig, chresp = multipath_fading_channel(tx_sig, pdp, n_dims)
    # rx_sig = repeat(sum(tx_sig, dims=2), 1, 2)
    rx_sig = awgn_channel!(rx_sig, N0)

    rx_frame = ofdm_demod(ofdm, rx_sig, tx_frame)
    # CFR = OFDM.Chest.convert_to_cfr(ofdm, chresp)
    # p = plot(abs.(CFR[:,1,1]))
    ave_pilots = pilot_symbol_average(ofdm, rx_frame)
    x = cir_estimate(ofdm, ave_pilots, N0=N0, method=:LMMSE)
    # CFRest = ave_pilots[:,1,1]
    # CFRest[ofdm.indices[:pilot]] ./ ofdm.pilot.symbols
end

# 変調パラメータ
N_FFT = 64
N_GI  = 16
N_DIMS=(2,2)
N_TX, N_RX = N_DIMS
BIT_LENGTH = 2000
M = 4
N_DATA = div(BIT_LENGTH, N_TX*Int(log2(M)))

# 変調パラメータ
modtype = Qam(M)
ofdm = Ofdm(n_fft=64, n_gi=16, n_dims=(2,2), pilot=(layout=:lte))
ofdmframe = OfdmFrame(ofdm, N_DATA)

# 通信路パラメータ
pdp = Uniform(interval=3, L=15) # PDP

tx_bits = bitrand(1024)
tx_data = mod(modtype, tx_bits)
tx_frame = deepcopy(ofdmframe)
