module LDPC

import FileIO
import Base: show, size

using LinearAlgebra
using SparseArrays: sparse, sparsevec, SparseVector, SparseMatrixCSC
using JLD2
using ..FEC: AbstractEncoder, AbstractDecoder


export LDPCEncoder,  
       LDPCDecoder

# source files
include("PEG.jl") # module
include("common.jl")
include("encode.jl")
include("decode.jl")



# 符号生成行列のオブジェクト保存先
const repos = joinpath(@__DIR__, "../../../../../data/ldpc_codes/")

# 生成行列, パリティ検査行列の読み込み
function load(src=repos; rows, cols, method="PEG", isregular=true, deg=0)
    prefix = ""
    prefix *= method

    if !isregular
        prefix = join([prefix, "irr", deg], "_")
    end

    rows = string(rows); cols = string(cols)
    fname = join([prefix, rows, cols], "_")

    fpath = repos * fname * ".jld2"

    local G, H
    try
        println(fname * " is now loading ...")
        if method=="PEG"
            rank = cols
        else
            rank = FileIO.load(fpath, "rank")
        end
        I, J = FileIO.load(fpath, "genmat")
        # G = row_list(Int32.(I),Int32.(J))
        G = sparse(Int32.(I), Int32.(J), ones(Int64,length(I)))
        I, J = FileIO.load(fpath, "paritymat")
        H = sparse(Int32.(I), Int32.(J), ones(Int64,length(I)))
    catch e
        error(e)
    end

    println("LDPC code has been successfully loaded.")
    return G, H
end


end # module
