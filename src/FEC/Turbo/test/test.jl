# Turbo.jl 動作確認
include(pwd()*"/Mymodule/Commun/Commun.jl")
Turbo = Commun.FEC.Turbo
using Random
n_bits = 1000
encoder, decoder = Turbo.TurboCoder(K=4, poly=[13 15 17], n_bits=n_bits)

txbits = bitrand(n_bits)
encbits = encoder(txbits)
txsyms = Commun.Mod.qammod(encbits, M=2)
rxsyms = Commun.awgn_channel(txsyms, N0=0.01)
lch = 4*rxsyms/0.01
