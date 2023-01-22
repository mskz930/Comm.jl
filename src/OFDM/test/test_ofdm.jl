# OFDM変復調のテスト
using Revise
using Comm.OFDM
using Random: bitrand

# パラメータ
pdp = ones(5); pdp ./= sum(pdp)

n_dims = (2,2)
n_bit = 4000
bps = 2
M = 2^bps
qam = Qam(M)
n_data = cld(n_bit, bps*n_dims[1])

pilotp = Pilot(:comb, dt=1, df=cld(64,8))
ofdmp = Ofdm(n_fft=64, n_gi=16, n_dims=n_dims, pilot=pilotp) # OFDMパラメータ
virtualmap!(ofdmp, n_data)
framelen = length(ofdmp.maplist)
frame = OfdmFrame(ofdmp, framelen, pilot_insert=true)

tx_frame = copy(frame)

# 送信データ
tx_bits  = bitrand(n_bit)
tx_syms  = mod(qam, tx_bits) # QAM変調
tx_sig = ofdm_mod(ofdmp, tx_syms, frame=tx_frame)

rx_sig = tx_sig
rx_frame = ofdm_demod(ofdmp, rx_sig, frame=tx_frame)
