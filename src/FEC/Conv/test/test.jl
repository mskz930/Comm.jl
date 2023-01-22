using Comm.Comutils
using Comm.FEC.Conv

# trellis = poly2trellis([133 171], fbpoly=133)
trellis = poly2trellis([13 17], fbpoly=13)
nbits = 100
decoder = APPDecoder(trellis, numbits=nbits)

tx_bits = randbit(nbits)
enc_bits = convenc(tx_bits, trellis)

rx_bits = mod(Qam(2), enc_bits)

Juno.@enter decoder(rx_bits)

trellis.outputs
