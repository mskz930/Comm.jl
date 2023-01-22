using Revise
using Commun.Comch
using Commun.Digimod
using Random

pdp = Uniform(space=1, τmax=0)

tx_bits = bitrand(1000)
tx_data = mod(Qam(4), tx_bits)

rx_data = similar(tx_data)
Juno.@enter multipath_fading_channel(tx_data, pdp, 1, fdTs=0.0, sps=1, tail=false)
