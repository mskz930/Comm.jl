using Comm
using .Modulator, .OFDM

# パラメタ
nfft = 64
ngi  = 16
ndim = (2,2)
pilot = Pilot("LTE", dt=1, numofpilots=16, style=:random)
ofdm = OFDM.Ofdm(nfft=nfft, ngi=ngi, ndim=ndim, pilot=pilot)
frame = OFDM.gen_frame(ofdm, 500)

#%%
# 送信データ生成
qam = Qam(4)
tx_bits = randbit(1000)
tx_data = mod(qam, tx_bits)
tx_frame = subcmap(ofdm, frame, tx_data)
tx_sig  = mod(ofdm, tx_frame)

# 通信路
pdp = [1, 0, 1] |> x -> x ./ sum(x)
rx_sig, chresp = Comch.multipath_fading(tx_sig, pdp, ndim)

# 受信機
rx_frame = OFDM.demod(ofdm, rx_sig)
CFR = OFDM.Chest.to_cfr(ofdm, chresp)

CFR_est  = OFDM.chest(ofdm, rx_frame, tx_frame)

rx_data = OFDM.equalize(ofdm, rx_frame, CFR)
rx_bits = demod(Qam(4), rx_data)
sum(tx_bits .!== rx_bits)
