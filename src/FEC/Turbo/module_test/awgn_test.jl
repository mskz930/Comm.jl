
include("../Turboc.jl")
include(joinpath(pwd(), "MyModule/Commun/Modulation/Mod.jl"))
include(joinpath(pwd(), "MyModule/Commun/Channel/Comch.jl"))
using Random: randperm, bitrand
using Plots

# グローバル定数
const TOTAL_BITS = 4e6 # 最大送信ビット数
const EPS = 1e4 # 誤りしきい値

# 送信ビット数
n_bits = 5000
# 符号化率
rate = 1//2
# ターボ符号パラメータ
turbo = Turbo.TurboCoder(K=4, poly=[13 15 17], fbpoly=13, n_bits=n_bits, iterations=6, rate=rate, method="max")
# 符号長
ncbits = Int(n_bits/rate) + turbo.n_tailbits
# 送信機インタリーバ
permind = randperm(ncbits)

# main
function main(ebnodBs)
    ebnodBs = collect(ebnodBs)
    ber_results = Float64[]
    for (i,ebnodB) in enumerate(ebnodBs)
        # Eb/N0修正
        m_ebnodB = ebnodB + 10*log10(turbo.rate)
        snrdB = m_ebnodB
        # 雑音電力密度
        N0 = 10^(-snrdB/10)
        # 初期化
        numbits=0; errors=0
        @time while numbits<TOTAL_BITS && errors<EPS
            global n_bits, turbo, nc_bits, permind
            # 送信ビット生成
            tx_bits = bitrand(n_bits)
            # ターボ符号化
            enc_bits = Turbo.encode(turbo, tx_bits)
            # インタリーブ
            enc_bits = enc_bits[permind]
            # QAM変調
            tx_symbols = Mod.qammod(enc_bits, M=2)
            # AWGN通信路
            rx_symbols = Comch.awgn!(tx_symbols, N0=N0/2)
            # チャネル尤度
            Lch = Mod.qamdemod(rx_symbols, N0, M=2)
            # デインタリーブ
            Lch[permind] = Lch
            # ターボ復号
            Lapp = Turbo.decode(turbo, Lch)
            # 受信LLR判定
            rx_bits = Lapp .> 0
            # 誤り数計算
            error = sum(tx_bits .!== rx_bits)
            # 送信ビット, 誤り数カウント
            errors += error; numbits += n_bits
        end
        # 平均ビット誤り率
        ber = errors/numbits
        print("Eb/N0: ", ebnodB, "[dB],    ")
        print("BER: ", ber)
        print(" (", numbits, "/",TOTAL_BITS, ")\n")
        if ber > 0
            push!(ber_results, ber)
        else
            println("プログラムを終了します...")
            ebnodBs = ebnodBs[1:i]
            break;
        end
    end
    Dict(:EBN0=>ebnodBs, :BER=>ber_results)
end

# test
function test(ebnodB)
    global n_bits, turbo, nc_bits, permind
    N0 = 10^(-ebnodB/10) # N0
    # 送信ビット生成
    tx_bits = bitrand(n_bits)
    # ターボ符号化
    enc_bits = Turbo.encode(turbo, tx_bits)
    # インタリーブ
    enc_bits = enc_bits[permind]
    # QAM変調
    tx_symbols = Mod.qammod(enc_bits, M=2)
    # AWGN通信路
    rx_symbols = Comch.awgn!(tx_symbols, N0=N0/2)
    # チャネル尤度
    Lch = Mod.qamdemod(rx_symbols, N0, M=2)
    # デインタリーブ
    Lch[permind] = Lch
    @time Lapp = Turbo.decode(turbo, Lch)
    # 受信LLR判定
    rx_bits = Lapp .> 0
    # 誤り数計算
    error = sum(tx_bits .!== rx_bits)
end
