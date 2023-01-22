using Revise, Random, Printf
using Commun
using Commun.Digimod

Random.seed!(1234)

nbit = 2048
M = 4
bps = Int(log2(M))
nfft, ngi, ndims = 64, 0, (2,3)
ndata = div(nbit, bps)

pdp = Exponent(ratio=20, L=8)
const qam = Qam{4}()
const ofdm = Ofdm(nfft, ngi, ndims)
const frame = OFDM.Frame(ofdm, ndata)

txframe = deepcopy(frame)

function test(SNR)
    global nbit, ndims, txframe, ofdm, qam, pdp
    N0 = 10^(-SNR/10)
    return transceiver(nbit, ndims, qam, ofdm, txframe, pdp, N0)
end

function transceiver(nbit, ndims, qam, ofdm, txframe, pdp, N0)
    txbits = Array(bitrand(nbit))
    txdata = Digimod.mod(qam, txbits)
    txsig  = OFDM.mod(ofdm, txdata, txframe)

    rxsig, chresp = multipath_fading_channel(txsig, pdp, ndims, N0)

    rxframe = OFDM.demod(ofdm, rxsig, txframe)
    cfr = OFDM.Chest.convert_to_cfr(ofdm, chresp)
    rxdata = OFDM.GSC.equalize(ofdm, rxframe, chresp, cfr, N0)
    rxbits = Digimod.demod(qam, rxdata)

    n_err_bit = sum(txbits .!= rxbits)
    return n_err_bit
end
