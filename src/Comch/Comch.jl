module Comch

using OffsetArrays
import StatsBase

export awgn!,
       awgn,
       rayfading,
       multipath_fading,
       multipath_fading!,
       gen_fading_process

include("PDP.jl")
include("transversal_filter.jl")
include("fading.jl")
include("channel.jl")

end
