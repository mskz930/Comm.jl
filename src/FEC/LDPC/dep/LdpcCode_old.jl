module LDPC
include("../LinearCode.jl")
include("./peg.jl") # module

using LinearAlgebra, SparseArrays
using .LinearCode

export PEG, ldpcdec, ldpcenc

abstract type Node end # ノード型定義


""" LDPC符号化 """
# 符号化オブジェクト
struct LdpcEncoder
    G::AbstractMatrix{Float64}
    rank::Integer
    function LdpcEncoder(G)
        if eltype(G)==Bool
            G = float.(G)
        end
        new(G)
    end
end

# 符号化関数
function ldpcenc(x::AbstractVector{}, G::AbstractMatrix{<:AbstractFloat})
    n = size(G,1) # 符号長
    k = size(G,2) # 情報ビット長
    x = float.(x) # 実数変換
    encoded_bit = [x; (view(G,k+1:n,:)*x).%2.0]
end

# 符号化関数(オブジェクトによる定義)
function (enc::LdpcEncoder)(x::AbstractVector)
    x = float.(x)
    n = size(enc.G,1) # 符号長
    k = size(enc.G,2) # 情報ビット長
    encoded_bits = Bool.([x; (view(enc.G,k+1:n,:)*x).%2.0])
end



""" LDPC復号化 """
# 復号器オブジェクト
struct LdpcDecoder
    H::AbstractMatrix{Bool}
    checknodes::Vector{<:Node}
    bitnodes::Vector{<:Node}
end
# 外部コンストラクタ
function LdpcDecoder(H::AbstractMatrix{Bool})
    bnlist, cnlist = get_indexlist(H)
    checknodes = [CheckNode(list) for list in bnlist]
    bitnodes = [BitNode(list) for list in cnlist]
    LdpcDecoder(H, checknodes, bitnodes)
end

# ノード番号のリストを行列から取得
function get_indexlist(H)
    bnlist = [LinearIndices(view(H,m,:))[view(H,m,:)] for m in axes(H,1)]
    cnlist = [LinearIndices(view(H,:,n))[view(H,:,n)] for n in axes(H,2)]
    return bnlist, cnlist
end


# チェックノード定義
struct CheckNode <: Node
    port::Dict{Int64,Float64}
    neighbor::Vector{Int64}
    function CheckNode(bnlist)
        port, neighbor = node_initialize(bnlist)
        new(port, neighbor)
    end
end

# ビットノード定義
struct BitNode <: Node
    port::Dict{Int64,Float64}
    neighbor::Vector{Int64}
    function BitNode(cnlist)
        port, neighbor = node_initialize(cnlist)
        new(port, neighbor)
    end
end

# Node初期化
function node_initialize(list::AbstractVector)
    port = Dict(i=>0.0 for i in list) # スパースベクトル初期化
    neighbor = list # リスト
    port, neighbor
end

# ldpcdec: LDPC復号
function ldpcdec(dec::LdpcDecoder, Lch::AbstractArray{<:AbstractFloat,1}, MAXITER=10, option=(:sum_product,:max); outtype=:soft)
    N = size(dec.H,2) # 符号長
    Lapp = zeros(Float64, N) # 最大事後尤度比
    decoded_bits = zeros(Bool, N) # 判定ビット配列
    # decoder_initialize!(dec.bitnodes, Lch)
    decoded_bits, Lapp, iter = log_sum_product(dec, decoded_bits, Lapp, Lch, MAXITER, option...)
    if outtype==:hard
        return decoded_bits
    elseif outtype==:soft
        return Lapp
    end
end

# オブジェクトによる関数定義
function (dec::LdpcDecoder)(Lch::AbstractArray{<:AbstractFloat,1}, MAXITER=10, option=(:sum_product,:max); outtype=:hard)
    N = size(dec.H, 2) # 符号長
    Lapp = zeros(Float64, N) # 最大事後尤度比
    decoded_bits = zeros(Bool, N) # 判定ビット配列
    # decoder_initialize!(dec.bitnodes, Lch)
    decoded_bits, Lapp, iter = log_sum_product(dec, decoded_bits, Lapp, Lch, MAXITER, option...)
    if outtype==:hard
        return decoded_bits
    elseif outtype==:soft
        return Lapp
    end
end

""" log_sum_product: sum-product/min-sum 復号 """
# log-sum-product
function log_sum_product(dec, decoded_bits, Lapp, Lch, MAXITER, method, cond)
    bitnodes = dec.bitnodes
    checknodes = dec.checknodes
    iter = 0
    while (iter < MAXITER)
        # bitnode to checknode
        message_passing!(Lch, bitnodes, checknodes,iter)
        # checknode to bitnode
        message_passing!(checknodes, bitnodes, method)
        if cond == :checksum
             decoded_bits = decision!(decoded_bits, Lapp)
             checksum = check_sum!(decoded_bits, dec.H)
            break
        end
        iter += 1
    end
    Lapp = lapp_calc!(Lapp, Lch, dec.bitnodes)
    decoded_bits = decision!(decoded_bits, Lapp)
    return decoded_bits, Lapp, iter
end

# message passing: bitnode to checknode
function message_passing!(Lch, from::Array{BitNode,1}, to::Array{CheckNode,1}, iter)
    Threads.@threads for n in eachindex(Lch) # ビットノード Threads.@threads
        neighbor = from[n].neighbor
        mes = Lch[n]
        if iter > 0
            for j in neighbor
                mes += from[n].port[j] # 外部メッセージ
            end
        end
        for i in neighbor
            to[i].port[n] = tanh((mes-from[n].port[i])/2.0)
        end
    end
end


# message passing: checknode to bitnode
function message_passing!(from::Array{CheckNode,1}, to::Array{BitNode,1}, method::Symbol=:sum_product)
    if method == :sum_product
        Threads.@threads for m in eachindex(from) # ビットノード
            neighbor = from[m].neighbor
            mes = 1.0 # 初期メッセージ
            for j in neighbor
                mes *= from[m].port[j]
            end
            for i in neighbor
                to[i].port[m] = 2.0*myatanh(mes./from[m].port[i])
            end
        end
    elseif method == :min_sum
        for m in eachindex(from) # ビットノード
            for m in from[m].neigbors
                alpha = 1.0; beta = 0.0
                for j in from[m].neighbors
                    if i !== j
                        val = from[m].port[j]
                        alpha *= sign(val)
                        beta  += phi(abs(val)/2)
                    end
                end
                to[i].port[m] = alpha*phi(beta)
            end
        end
    end
end

function myatanh(x)
    if x < 1.0
        y = atanh(x)
    else
        y = 19.07
    end
    return y
end

# phi function
function phi(x)
    if 0 < x
        y = log(coth(x/2))
    elseif x == 0
        y = convert(eltype(x),707.58)
    else
        error("xは正の実数でなければいけません。")
    end
end

# lapp_calc: 現在の事後比を集計する
function lapp_calc!(Lapp::AbstractVector, Lch,  bitnodes::Array{BitNode,1})
    for n in eachindex(bitnodes)
        Lapp[n] = Lch[n]
        for i in bitnodes[n].neighbor
            Lapp[n] += tanh(bitnodes[n].port[i]/2.0)
        end
    end
    return Lapp
end

# decision!: 判定関数 """
function decision!(x̂, Lapp)
    @inbounds for i in eachindex(Lapp)
        x̂[i] = Lapp[i] .< 0
    end
    return x̂
end

# checksum: チェックサムの計算 """
function check_sum!(Lapp, H)
    x̂ = Lapp .< 0 # 判定ビット
    checksum = false
    for j in axes(H,2)
        for i in axes(H,1)
            if H[i,j]
                checksum ⊻= H[i,j] & x̂[i]
            end
        end
    end
    return checksum
end


end # end module
