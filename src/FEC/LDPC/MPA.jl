module MPA

import Base.show
using SparseArrays
using LinearAlgebra: rank


export Decoder, LogSPA, MinSum


struct LogSPA end
struct MinSum end

"""
mutable struct Node
    val::Float64
    port::Dict{Int32,Float64}
    neighbors::Dict{Int32,Node}

ノード
"""
mutable struct Node
    val::Float64
    edges::Dict{Int32, Float64}
    neighbors::Dict{Int32,Node}
    Node() = new(0.0, Dict(), Dict())
end

show(io::IO, n::Node) = print(io, """$(typeof(n))((val->$(n.val), edges->$(keys(n.edges)))""")

show(io::IO, m::MIME"text/plain", n::Node) = print(io, """
   $(typeof(n)):
       val  -> $(n.val)
       edges -> $(keys(n.edges))
""")


"""
    Decoder(H)
メッセージパッシング復号器オブジェクト
"""
mutable struct Decoder
    # H::SparseMatrixSCS
    input_size::Int32
    output_size::Int32
    checknodes::Vector{Node}
    bitnodes::Vector{Node}
    method::Union{LogSPA, MinSum}
    max_iter::Int16

    function Decoder(H; method=:LogSPA, max_iter=30)
        m, n = size(H)
        checknodes = [Node() for _ in 1:m]
        bitnodes   = [Node() for _ in 1:n]
        bnlist, cnlist = makelist(H)
        _link!(checknodes, cnlist, bitnodes, bnlist)
        method = eval(Meta.parse(string(method)*"()"))
        new(k, n, checknodes, bitnodes, method, max_iter)
    end
end

function show(io::IO, d::Decoder)
    print(io, """(m,n)=($(length(d.checknodes)), $(length(d.bitnodes)))""")
end
function show(io::IO, m::MIME"text/plain", d::Decoder)
    print(io, """
    $(typeof(d))(
        (m,n) = ($(length(d.checknodes)), $(length(d.bitnodes)))
    )
    """)
end

# ノードリスト集合をパリティ検査行列か作成する
function makelist(H)
    bnlist = [findall(!iszero, H[m,:]) for m in axes(H,1)]
    cnlist = [findall(!iszero, H[:,n]) for n in axes(H,2)]
    return bnlist, cnlist
end

# 初期化して隣接ノードを保存する
function _link!(checknodes, cnlist, bitnodes, bnlist)
    # ノード同士を結合
    for (m,cn) in enumerate(checknodes)
        for n in bnlist[m]
            cn.edges[n] = 0.0
            cn.neighbors[n] = bitnodes[n]
        end
    end
    for (n,bn) in enumerate(bitnodes)
        for m in cnlist[n]
            bn.edges[m] = 0.0
            bn.neighbors[m] = checknodes[m]
        end
    end
end

# ポートの初期化
function _init!(d::Decoder)
    for m in eachindex(d.checknodes)
        for n in d.bnlist
            d.checknodes[m].edges[n] = 0.0
        end
    end
    for n in eachindex(d.bitnodes)
        for m in d.cnlist
            d.bitnodes[n].edges[m] = 0.0
        end
    end
end

"""
    (d::Decoder)(Lin)

LDPC復号
    args:
      Lin:
    returns:
      Lout:
"""
function (d::Decoder)(Lin, target=:output)
    len = target == :output ? d.output_size : d.input_size
    Lout = Array{Float64}(undef, d.output_size)
    checknodes, bitnodes, method, max_iter = d.checknodes, d.bitnodes, d.method, d.max_iter
    decode!(Lout, Lin, checknodes, bitnodes; method=method, max_iter=max_iter)
end

"""LDPC符号を復号する"""
function decode!(Lout, Lin, checknodes, bitnodes; method, max_iter)
    _init!(bitnodes, Lin) # 初期化

    n_iter = 0
    while n_iter < max_iter
        _to_checknode!(bitnodes, n_iter)
        _to_bitnode!(checknodes)
        _message_update!(bitnodes)
        cks = _checksum(checknodes)
        cks == 0 && break
        @show n_iter
        n_iter += 1
    end
    _copyto!(Lout, bitnodes)
    Lout
end


# 入力対数尤度を初期値としてセットする
function _init!(bitnodes, Lin)
    for n in eachindex(bitnodes)
        bitnodes[n].val = Lin[n]
    end
end

# チェックノード => ビットノード
function _to_checknode!(bitnodes, n_iter)
    if n_iter == 0
        for (n,bn) in enumerate(bitnodes)
            for (m,cn) in bn.neighbors
                cn.edges[n] = tanh(bn.val/2)
            end
        end
    else
        for (n,bn) in enumerate(bitnodes)
            for (m,cn) in bn.neighbors
                mes = tanh((bn.val - bn.edges[m])/2) # メッセージ
                mes = ifelse(mes>1e-20, mes, 1e-20)
                cn.edges[n] = mes
            end
        end
    end
end

# チェックノード => ビットノード
function _to_bitnode!(checknodes)
    for (m,cn) in enumerate(checknodes)
        mes = reduce(*, values(cn.edges))
        for (n,bn) in cn.neighbors
            bn.edges[m] = 2*atanh(mes / cn.edges[n])
        end
    end
end

# ビットノードの値を更新
function _message_update!(bitnodes)
    for bn in bitnodes
        bn.val = val = reduce(+, values(bn.edges))
    end
end

# チェックサムの計算
function _checksum(checknodes)
    cks = Bool(0)
    for cn in checknodes
        for bn in values(cn.neighbors)
            cks ⊻= bn.val < 0
        end
        cks > 0 && return cks
    end
    return cks
end

# 結果を配列に保存する
function _copyto!(Lout, bitnodes)
    for i in eachindex(Lout)
        Lout[i] = bitnodes[i].val
    end
    Lout
end


end # module
