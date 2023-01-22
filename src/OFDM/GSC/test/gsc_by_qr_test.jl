#%%
using Revise, Random, Plots
using Commun.OFDM
using Commun.OFDM.GSC
using Commun.OFDM.Chest: convert_to_cfr
using Commun.Comch
using Commun.Digimod

# オブジェクト
n_dims = (2,3)
modtype = Qam(4)
ofdmp = Ofdm(n_fft=64, n_gi=0, n_dims=n_dims, pilot=PilotType())
ofdmframe = OfdmFrame(ofdmp, framelen=10)

# 通信路パラメータ
pdp = zeros(7); pdp[[1,3,5]] .= 1.0; pdp ./= sum(pdp)
n_bit = 1000
SNR = 30
N0  = 10^(-SNR/10)
tx_frame = deepcopy(ofdmframe)

# 送信側
tx_bits = bitrand(n_bit)
tx_data = mod(modtype, tx_bits)
tx_sig  = ofdm_mod(ofdmp, tx_frame, tx_data)

# 通信路
rx_sig, chresp = multipath_fading_channel(tx_sig, pdp, n_dims)
#%%
rx_data = GSC.equalize_by_qr(ofdmp, rx_sig, chresp, N0, method=:SIC, slice=x->slice(modtype, x))

datainds = get_data_indices(tx_frame)
demod_data = rx_data[:,datainds]


rx_bits = demod(modtype, demod_data) |> vec

sum(tx_bits .!= rx_bits)
