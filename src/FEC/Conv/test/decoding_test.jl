""" This is a script file """
# 畳み込み符号 エンコード/デコード テスト
include("../conv.jl")
using Random
Random.seed!(1)
nbits = 1000 # 送信ビット数
ebnodB = 5 # Eb/N0[dB]
N0 = 10^(-ebnodB/10) # 雑音電力密度
tx_data = rand(nbits) .> 0.5 # バイナリデータ生成
trellis = conv.poly2trellis(3, [7 5]) # トレリス
encoded_data = conv.convenc(trellis, tx_data, isterm=true) # 符号化ビット列
tx_symbol = 2 .* encoded_data .- 1 # 送信シンボル(BPSK)
rx_symbol = tx_symbol .+ (N0/2)*randn(size(encoded_data)) # 受信シンボル
y = rx_symbol .> 0
y = reshape(y, trellis.inout[2], :)
conv.vithard(y, trellis, isterm=true)
