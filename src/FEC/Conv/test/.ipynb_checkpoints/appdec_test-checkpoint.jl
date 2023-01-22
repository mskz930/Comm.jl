""" BCJR　MAP復号テスト """
# include("./Commun/function/error_correcting/conv/test/test.jl")
# include("./Commun/function/error_correcting/conv/conv.jl")
include("../conv.jl")
using BenchmarkTools
using Random
Random.seed!(1);

trellis = conv.poly2trellis(3, [7 5], 7); # トレリス
mapdec = conv.MAPDecoder(trellis); # MAP復号器パラメータ
nbits = 1000; # 送信ビット数
ebnodB = 5
snrdB = ebnodB - log10(2) # 等価SNR
N0 = 10^(-snrdB/10); # 雑音電力密度
tx_data = rand(nbits) .> 0.5 # バイナリデータ生成
encoded_data = conv.convenc(trellis, tx_data, isterm=false);
tx_symbol = 2. .* encoded_data .- 1.; # 送信シンボル
rx_symbol = tx_symbol .+ (N0/2)*randn(size(encoded_data)); # 受信シンボル
llrbits = 4/N0 .* rx_symbol; # LLR bits
Lx = reshape(llrbits, trellis.inout[2],:)
La = zeros(eltype(Lx), trellis.inout[1], size(Lx,2))
Lapp = mapdec(Lx, La, "input")
decoded_data = Lapp[:] .> 0
sum(tx_data .!== decoded_data) # ビット誤り数
