module DigitalModulation

using Base: Bool
import Base: mod
using UnPack
using DataStructures: OrderedDict
using Binary: grayenc

include("common.jl")
include("pam.jl")
include("qam.jl")
include("psk.jl")

include("slicer.jl")
# include("theory.jl")

ModType = Union{Qam,Pam}

export Pam, Qam,
  SoftSymbol,
  ModType,
  slice,
  mod,
  mod!,
  smap,
  sdemap,
  demod,
  demod!,
  scaling, rescaling,
  normfactor,
  get_refs,
  make_mapping_table,
  make_refs
# make_bsets # methods


end # module
