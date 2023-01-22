using Comm
using .OFDM

function test(; ndim, method, args=(θ=0.3, n_list=1))
    qam = Qam(4)
    pilot = Pilot("comb", t0=1, dt=2)
    ofdm = Ofdm(nfft=64, ngi=16, ndim=(2,2), pilot=pilot)
    frame = gen_frame(ofdm,1400)
    N0 = 0.01

    tx_data = rand(qam, 1400)
    tx_frame = subcalloc(ofdm, frame, tx_data)
    tx_sig  = ofdm_mod(ofdm, tx_frame)

    pdp = Comch.PDP.exponent(L=8, ratio=20)
    rx_sig, chresp = multipath_fading(tx_sig, pdp, (2,2), N0)
    rx_frame = ofdm_demod(ofdm, rx_sig)
    CFR = Chest.to_cfr(ofdm, chresp)
    CFR_est = Chest.chest(ofdm, rx_frame, tx_frame, interp=false, ave=true)
    CIR_est = Chest.cirest(ofdm, CFR_est, N0, L=ofdm.ngi, method=method, args=args)
end
