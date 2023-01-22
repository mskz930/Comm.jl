# OFDM.jl
module OFDM

# dependencies
import Base: mod
using FFTW: fft!, fft, ifft!, ifft, plan_fft, plan_ifft
using Random: randperm, shuffle!, shuffle
using Statistics: mean
using StatsBase: sample
using UnPack
using LinearAlgebra: Diagonal
using Myutils: shift!, zero_padding

using My
using ..DigitalModulation
using ..SignalDetection
using ..Comutils: dftmat

export Pilot,
       make_time_table,
       make_state_table,
       get_time_idxs,
       getdim,
       Comb, Block, LTE, Scatter, Layout

export Frame,
       Ofdm,
       OfdmMod,
       CP,
       ZP,
       genframe,
       ofdm_mod,
       ofdm_demod,
       subcmap,
       subcmap!,   #
       mapping,
       demapping,  # データシンボルの抽出
       dftsmap,
       dftsdemap,
       serial_to_parallel,
       parallel_to_serial,
       get_idxs,
       get_idxs_list,
       extract_data,
       equalize,            # チャネル等化(信号検出)
       Chest,
       CPR


# source file
include("utils.jl")         # 共通関数
include("preamble.jl")      # プリアンブル信号
include("pilot.jl")         # パイロット
include("dataframe.jl")         # 信号フレーム
include("OfdmMod.jl")      # パラメタ型
include("tools.jl")         # OFDM処理用自作関数
include("mapping.jl")       # サブキャリアマッピング
include("modulation.jl")    # OFDM変調
include("equalize.jl")      # チャネル等化
# include("chest.jl")         # チャネル推定
include("Utils/Utils.jl")   # Utility

# module file
include("Chest/Chest.jl")   # パイロットチャネル推定
include("CPR/CPR.jl")       # CPR equalizer





end # module
