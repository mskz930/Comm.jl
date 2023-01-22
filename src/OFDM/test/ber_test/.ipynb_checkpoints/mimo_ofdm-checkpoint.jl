module Test
include("../../ofdm.jl")
include("../../../../Channel/Channel.jl")
include("../../../../Utility/Conversions.jl")

using  LinearAlgebra, FFTW
using .Conversions, .Ofdm, .Channel

# グローバルパラメータ
nbits = 1000 # 送信ビット数
nfft=64; cpsize=0; # FFT数, CP長
npi=12; ngc=6; # パイロットキャリア, ガードキャリア
ndim=(2,2) # 送受信アンテナ数
params = OfdmParams(nfft,cpsize, npi, ngc, ndim, pilot_type=:comb, pilot_space=4, pilot_interval=2) # OFDMパラメータ
indices = get_indices(params) # インデックス
frame = gen_frame(params, indices, 30); # フレームテンプレート
nbits = 1000
M = 4; bps = Int(log2(M)) # 送信ビット数
qam = Ofdm.QamMod(M) # 変調器
pdp = [1.0]; # PDP
pdp ./= sum(pdp) # 正規化
pdp = Channel.exponent(5,20,space=2);

function main(snrdBs)
    ser_results = Float64[]
    for snrdB in snrdBs
        MAXITER = 0
        N0 = 10^(-snrdB/10) # nvar
        errors = 0; numsyms = 0
        while (errors < 4e3) && (numsyms < 4e5)
            tx_ofdmsig, tx_bit = transmitter(nbits) # 送信機
            rx_ofdmsig, cir = channel(tx_ofdmsig, pdp, N0) #　通信路
            rx_bit = receiver(rx_ofdmsig, cir, MAXITER) # 受信機
            errors += symerror(tx_bit, rx_bit, bps)
            numsyms += nbits/bps
        end
        ser = errors/numsyms
        push!(ser_results,ser)
        @show snrdB, ser
    end
    return snrdBs, ser_results
end

function test(snrdB)
    N0 = 10^(-snrdB/10) # 雑音電力密度
    MAXITER = 1
    tx_ofdmsig, tx_bit = transmitter(nbits) # 送信機
    rx_ofdmsig, cir = channel(tx_ofdmsig, pdp, N0) # 通信路
    rx_bit= receiver(rx_ofdmsig, cir, MAXITER) # 受信機
    errors = biterror(tx_bit, rx_bit)
end

# 送信機
function transmitter(nbits, frame=frame, params=params, indices=indices)
    tx_bit = rand(nbits) .> 0.5 # 送信ビット
    tx_symbol = Ofdm.qammod(qam, tx_bit) # 送信シンボル
    tx_frame = Ofdm.data_mapping(params, frame, tx_symbol, indices); # 送信フレーム
    tx_ofdmsig = Ofdm.ofdmmod(params, frame, tx_symbol, indices); # 送信OFDM信号
    return tx_ofdmsig, tx_bit
end

# 通信路
function channel(tx_ofdmsig, pdp, N0, ndim=ndim)
    rx_ofdmsig, cir = Channel.multipath_fading(tx_ofdmsig, ndim, pdp, 0.0, istail=false);
    Channel.awgn!(rx_ofdmsig, N0)
    return rx_ofdmsig, cir
end

# 受信機
function receiver(rx_ofdmsig, cir, MAXITER, temp_frame=copy(frame), params=params, indices=indices)
    frameinfo = indices[:frame]
    rx_frame = ofdmdemod(params, rx_ofdmsig); # 受信フレーム
    cfr_temp = get_cfr(cir, 64, 0); # 周波数応答
    cfr = repeat(cfr_temp, 1, size(rx_frame,2), 1, 1) #
    nblock = size(rx_frame, 2) # OFDMシンボル数
    buffer = eltype(rx_ofdmsig)[]; sizehint!(buffer, nblock*params.Nt*params.data_carriers[1]) # バッファ
    rind = [LinearIndices(view(frameinfo,:,j,:))[view(frameinfo,:,j,:) .!== 0] for j in 1:nblock] # ヌル以外のインデックス
    dind = [LinearIndices(view(frameinfo,:,j,:))[view(frameinfo,:,j,:) .!== 0] for j in 1:nblock] # データインデックス
    DFTmtx = fft!(Array{ComplexF64}(I,params.nfft, params.nfft), 1)
    iter = 0
    # CPR部分
    while iter <= MAXITER
        eq_frame = copy(rx_frame) # 等化用フレーム
        # チャネル推定
        if iter > 0
            #=
            # 前判定シンボルを用いてパイロットシンボルのCPR
            for j in size(rx_frame, 2)
                # チャネル推定値補正
                current_symbol = view(temp_frame, :, j, :)
                cpr_correct!(cfr_temp, current_symbol)
                cfr = repeat(cfr_temp,1,size(rx_frame,2),1,1)
            end
            # チャネル推定
            # cir_est = cir_estimation() # CIRの推定
            =#
        end
        # CPR
        for j in 1:nblock
            if iter==0 && j==1
                previous_symbol = nothing
                current_symbol = nothing
            elseif iter==0 && j>1
                previous_symbol = view(temp_frame,:,j-1,:)
                current_symbol = nothing
            elseif iter>0 && j==1
                previous_symbol = nothing
                current_symbol = view(temp_frame,:,j,:)
            elseif iter>0 && j>1
                previous_symbol = view(temp_frame,:,j-1,:)
                current_symbol = view(temp_frame,:,j,:)
            end
            previous_symbol = nothing
            CPR.cpr!(view(eq_frame,:,j,:), cir, previous_symbol, current_symbol, DFTmtx, params.nfft, params.cpsize, rind[j], dind[j]) # 干渉除去
            eq_symbol = equalize(view(eq_frame,:,j,:), view(cfr,:,j,:,:), view(frameinfo,:,j,:)) # 信号検出
            temp_symbol = Ofdm.slicer(qam, eq_symbol) # スライス
            data_remapping!(params, view(temp_frame,:,j,:), temp_symbol, view(frameinfo,:,j,:)) # データマップ
            if iter == MAXITER
                buffer = append!(buffer, eq_symbol) # バッファに等化シンボルを溜める
            end
        end
        iter += 1
    end
    rx_bit = Ofdm.qamdemod(qam, buffer)
    return rx_bit
end

end
