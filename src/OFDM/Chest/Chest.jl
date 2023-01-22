module Chest

using Statistics: mean, mean!
using FFTW: fft, fft!, ifft, ifft!
using LinearAlgebra: pinv, inv, I, Diagonal, norm, diag
using UnPack

using CompressedSensing
# import .Greedy, .BCS.BMP, .BCS.SBL

using ..OFDM
using ....Comutils: dftmat
using Myutils: squeeze

# structs
struct Time end
struct Freq end


# source files
include("utils.jl")             # 汎用関数
include("average.jl")           # 平均化
include("interpolate.jl")       # 補間用関数
include("cir.jl")               # インパルス応答推定
include("cfr.jl")               # 周波数応答推定

# include("channe_estimate.jl")

cfrest(ofdm::OfdmMod, rx_frame, tx_frame) = channel_estimate(rx_frame, tx_frame; ofdm.chest...)

# チャネル推定(配列がある場合)
function cfrest(rx_frame, tx_frame; interp=true, exterp=true, average=false, domain=[:freq])
    pilot = ofdm.pilot
    isnothing(pilot) && error("Pilotが設定されていません。")

    cfr = zeros(eltype(tx_frame), size(rx_frame)..., size(tx_frame,3))
    
    if divide
        cfr = _divide!(cfr, rx_frame, tx_frame, ofdm.frame.idxs)
    else
        cfr = extract_pilot!(cfr, rx_frame, ofdm.frame.idxs)
    end
end

# パイロット除算
function _divide!(cfr, rx_frame, tx_frame, idxs)
    for j in axes(cfr,2)
        for i in axes(cfr,1)
            p = idxs[i,j]
            p <= 0 && continue
            for q in axes(rx_frame,3)
                cfr[i,j,q,p] = rx_frame[i,j,q] / tx_frame[i,j,p]
            end
        end
    end
    cfr
end



# L: maximum channel length
function noise_reduction(CFR, L)
    N = size(CFR,1)
    shape = size(CFR)
    Fl = dftmat(N)[:, 1:L] ./ sqrt(N) # (N x L) DFT matrix
    CFR = reshape(CFR, N, :)
    CFR = Fl * ((Fl'*Fl) \ (Fl'CFR))
    CFR = reshape(CFR, shape...)
    CFR
end


# 観測行列を作成
function get_dft_matrix(ofdm::Ofdm, L=ofdm.n_gi)
    n_fft, pidxs = ofdm.n_fft, ofdm.pilot.idxs
    F = dftmat(n_fft)[pidxs,1:L+1]
    F
end

function get_measurement_matrix(ofdm::Ofdm, L=ofdm.n_gi)
    @unpack n_fft, pilot = ofdm
    pidxs = pilot.idxs
    X = Diagonal(pilot.symbols)
    F = dftmat(n_fft)[pidxs,1:L+1]
    A = X * F
    A
end



end # module
