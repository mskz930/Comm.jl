""" BCJR　MAP復号テスト """
# include("./Commun/function/error_correcting/conv/test/mapdec_test.jl")
# %%
include("../conv.jl")
function ber_calc(ebnodBs)
    nbits = 1000; # 送信ビット数
    trellis = conv.poly2trellis(3, [7 5], 7); # トレリス
    mapdec = conv.MAPDecoder(trellis); # MAP復号器パラメータ
    ber_results = zeros(Float64, length(ebnodBs))
    for (n, ebnodB) in enumerate(ebnodBs)
        println("ebnodB:", ebnodB)
        snrdB = ebnodB - log10(2)
        N0 = 10^(-snrdB/10); # 雑音電力密度
        trans_bits = error_bits = 0
        while (trans_bits < 4e7)
            tx_data = rand(nbits) .> 0.5 # バイナリデータ生成
            encoded_data = conv.convenc(trellis, tx_data, isterm=false); # 符号化列
            tx_symbol = 2 .* encoded_data .- 1; # 送信シンボル
            rx_symbol = tx_symbol .+ (N0)*randn(size(encoded_data)); # 受信シンボル
            llrbits = 4/N0 .* rx_symbol; # LLR bits
            Lx = reshape(llrbits, trellis.inout[2],:) # 入力LLR
            La = zeros(eltype(Lx), trellis.inout[1], size(Lx,2)) # 事前LLR
            Lapp = mapdec(Lx, La)
            decoded_data = Lapp[:] .> 0
            error_bits += sum(tx_data .!== decoded_data) # ビット誤り数
            trans_bits += nbits # 送信ビット総数
        end
        ber_results[n] = error_bits / trans_bits # ビット誤り率
    end
    return ber_results
end

using Random, Plots
Random.seed!(1);
ebnodBs = 0; # EbN0dBs
@time ber_results = ber_calc(ebnodBs)
# %%
plot(ebnodBs, ber_results, yaxis = ((1e-7, 1e0), :log))
