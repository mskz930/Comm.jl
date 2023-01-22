module Test
include("../../Ofdm.jl")
include("../../../../Channel/Channel.jl")
include("cpr_receiver.jl")
using .Ofdm, .Channel
using FFTW, LinearAlgebra
nbits = 300; # 情報ビット数
M = 2; # 変調サイズ
qam = QamMod(M) # Qam変調オブジェクト
iterations = 1
nfft=64; cpsize=0; npi=0; ngc=0; ndim=(1,1) # OFDMパラメータ定義
params = OfdmParams(nfft,cpsize,npi,ngc,ndim) # OFDMパラメータ
indices = get_indices(params) # index集合
data_frame = gen_frame(params, indices, 500); # 送信データフレーム
pdp = zeros(5) # pdp
pdp[[1,3,end]] .= 1.0; pdp ./= sum(pdp) # 3-pathmodel

function test()
    tx_bit = rand(nbits) .> 0.5
    tx_symbol = qammod(qam, tx_bit) # 送信シンボル
    tx_frame = data_mapping(params, data_frame, tx_symbol, indices[:frame])
    tx_ofdmsig = ofdmmod(params, data_frame, tx_symbol, indices) # 送信OFDM信号
    rx_ofdmsig, cir = multipath_fading(tx_ofdmsig, pdp, ndim; fdTs=0.0, sps=params.Ts) # チャネル応答
    rx_frame = ofdmdemod(params, rx_ofdmsig) # 受信フレーム
    # ave_cir = Ofdm.get_avecir(cir, params.nfft, params.cpsize, size(rx_frame,2)) # 平均的なチャネル応答取得
    cfr = fft(vcat(cir, zeros(params.nfft-size(cir,1), size(cir)[2:end]...)),1) # 周波数応答
    eq_frame, nvar = hard_cpr(rx_frame, tx_frame, params, indices, iterations; N0=0.0, cir=cir, cfr=cfr, mod=qam) # 受信シンボル
end

end
