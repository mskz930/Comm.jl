using Comm
using .OFDM, .Modulator
using .Comutils: squared_error

#%%
SNR = 10 # dB
N0 = 10^(-SNR/10)
nbits = 4000
nfft, ngi = 128, 16
qam = Qam(4)
pilot = Pilot("LTE", t0=1, dt=6, df=7)
ofdm = Ofdm(nfft=nfft, ngi=ngi, ndim=(2,2), ngc=8, pilot=pilot)
sigframe = genframe(ofdm, div(nbits,qam.bps))


tx_data = mod(qam, randbit(nbits))
tx_frame = subcmap(ofdm, sigframe, tx_data)
tx_sig  = ofdm_mod(ofdm, tx_frame)

pdp = Comch.PDP.exponent(L=8, ratio=20)
rx_sig, chresp = multipath_fading(tx_sig, pdp, (2,2), N0)
rx_frame = ofdm_demod(ofdm, rx_sig)

#%%
CFR = Chest.to_cfr(ofdm, chresp)
CFRest = Chest.channel_estimate(ofdm, rx_frame, tx_frame, interp=:linear, average=:time, exterp=:flat)
CIR, _ = Chest.cirest(ofdm, CFRest, method=:interp_LS)

squared_error(chresp, CIR)


sticks(abs.(reshape(chresp, :, 4)), marker=true, layout=(2,2))
sticks!(abs.(reshape(CIR, :, 4)), marker=true, layout=(2,2))

plot(abs.(CFR[:,1,2]), marker=true)
plot!(abs.(CFRest[:,1,2]), marker=true)
plot!(abs.(Chest.to_cfr(ofdm, CIR[:,1,2])), marker=true)
