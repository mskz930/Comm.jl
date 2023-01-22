module Conv

using Binary
using SparseArrays
using UnPack
using DataStructures
using ..FEC: AbstractEncoder, AbstractDecoder

export convenc,
       convenc!,
       puncture,
       depuncture,
       get_punc_mtx,
       Trellis,
       poly2trellis,
       code_rate,
       num_of_outputs,
       make_transmat,
       BCJRDecoder,
       APPDecoder

# source files
include("utils.jl")
include("trellis.jl") 
include("encode.jl") 
include("decode.jl")
include("viterbi.jl")

end
