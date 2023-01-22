module FEC

abstract type AbstractEncoder end
abstract type AbstractDecoder end

# module files
include("LinearCodes.jl")  # 線形符号
include("Conv/Conv.jl")   # 畳み込み符号
include("Turbo/Turbo.jl") # ターボ符号
include("LDPC/LDPC.jl")   # LDPC符号

export Conv, LDPC
export AbstractEncoder, AbstractDecoder


end
