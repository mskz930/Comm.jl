# appdecoder_test.jl
using Comm
using Comm.FEC.Conv

# オブジェクト生成
trellis = poly2trellis([171 133], fbpoly=171)
appdec  = APPDecoder(trellis, bitlen=100)

tx_bits  = 
enc_bits = convenc(tx_bits, trellis)
