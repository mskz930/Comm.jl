module Comm

include("Utils/Utils.jl")
include("Comutils/Comutils.jl")
include("FEC/FEC.jl")
include("DigitalModulation/DigitalModulation.jl")
include("Comch/Comch.jl")
include("SignalDetection/SignalDetection.jl")
include("OFDM/OFDM.jl")

Digimod = DigitalModulation # alias
SigDet = SignalDetection

using .Comutils
using .Comch
using .DigitalModulation
using .OFDM
using .FEC

export Utils,
  Comutils,
  FEC,
  DigitalModulation, 
  Digimod,
  Comch,
  SignalDetection, 
  SigDet,
  OFDM

export randbit, awgn, multipath_fading, crandn

export QamMod, Qam


end # module
