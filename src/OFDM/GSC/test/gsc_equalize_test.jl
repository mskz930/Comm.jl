# GSC.equalizeのtest
using Random
using Plots
using Commun

nfft, ngi = 64, 16
ndim = (1,2)
nbit = 2048
M = 4; bps = ndigits(M-1, base=2)
ndata = div(nbit, bps * ndim[1])
const qam = Digimod.Qam{M}()
const ofdm = OFDM.Ofdm(nfft, ngi, ndim, pilot=:none)
const ofdmframe = OFDM.OfdmFrame(ofdm, ndata)
tx_frame = deepcopy(ofdmframe)
pdp = Comch.Exponent(ratio=10, L=12)

# 送受信テスト
function test(SNR)
    global nbit, ndim, qam, ofdm, pdp, tx_frame
    N0 = 10^(-SNR/10)

    # 送信側
    tx_bits  = bitrand(nbit)
    tx_data = Digimod.mod(qam, tx_bits)
    tx_sig  = OFDM.mod(ofdm, tx_data, tx_frame)

    # 通信路
    rx_sig, chresp = Comch.multipath_fading_channel(tx_sig, pdp, ndim)
    Comch.awgn_channel!(rx_sig, N0)

    # GSC受信機
    rx_frame = OFDM.demod(ofdm, rx_sig, tx_frame)
    cfr = OFDM.Chest.convert_to_cfr(ofdm, chresp)
    if ofdm.n_gi < size(chresp,1)-1
        rx_data  = OFDM.GSC.equalize(ofdm, rx_frame, chresp, cfr, N0)
    else
        rx_data = OFDM.equalize(ofdm, rx_frame, cfr, method=:MRC)
    end
    p = scatter(rx_data)
    display(p)
end
