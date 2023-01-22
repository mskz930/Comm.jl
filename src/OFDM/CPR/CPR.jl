module CPR

using FFTW, LinearAlgebra
using UnPack
using ..OFDM
using ..DigitalModulation: SoftSymbol
using Base.Iterators: flatten

include("utils.jl")              # 汎用関数
include("coef_matrix.jl")        # 係数行列生成
include("cpr_in_time.jl")        # Weigted CPR
include("cpr_in_freq.jl")        #
# include("soft_cpr_in_freq.jl")   # Soft CPR
# include("kalman_equalizer.jl")   #
# include("cp_reconstruct.jl") # cpr関数
# include("cpr_receiver.jl")




end # end module
