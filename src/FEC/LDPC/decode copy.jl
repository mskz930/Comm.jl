using UnPack

abstract type Algorithm end
struct LogSumProduct <: Algorithm end
struct MinSum <: Algorithm end


struct Table
  rows::Vector{SparseVector{Float64,Int64}}
  cols::Vector{SparseVector{Float64,Int64}}
end

struct Decoder <: AbstractDecoder
    H::SparseMatrixCSC
    bnlist::Vector{Vector{Int64}}                # 周辺ビットノードリスト
    cnlist::Vector{Vector{Int64}}                # 周辺チェックノードリスト
    rows::Vector{SparseVector{Float64,Int64}}    # メッセージ交換用擬似行列-行
    cols::Vector{SparseVector{Float64,Int64}}    # メッセージ交換用擬似行列-列
    params::Dict{Symbol, Int}                    # 復号用パラメータ
    method::Symbol                               # 復号メソッド

    function Decoder(H; max_iter=30, method=:log_sum_product)
        bnlist, cnlist = _get_indexlist(H)
        rows    = [sparsevec(i, zeros(Float64, length(i))) for i in bnlist]
        cols = [sparsevec(i, zeros(Float64, length(i))) for i in cnlist]
        params = Dict(:max_iter=>max_iter, :cks=>0)
        new(H, bnlist, cnlist, rows, cols, params, method)
    end
end

Base.show(io::IO, d::Decoder) = print(io, "$(typeof(d)): $(size(d.H)) LDPC Decoder")
Base.show(io::IO, m::MIME"text/plain", d::Decoder) = print(io, """$(typeof(d)): """)

# ノード番号のリストを行列から取得
function _get_indexlist(H)
    bnlist = @views [findall(!iszero, H[m,:]) for m in axes(H, 1)] # bit node list
    cnlist = @views [findall(!iszero, H[:,n]) for n in axes(H, 2)] # check node list
    bnlist, cnlist
end

"""
LDPC decoding by sum-product-algorithm
"""
function (dec::Decoder)(inp; flip=true, target=:output)
  out = decode(d, inp; flip=flip, dec.params...) # APP-LLRs
  if target == :output
      return out
  else
      return @views out[1:d.n_inputs]
  end
end



"""
  log_sum_product(Lpos, Lpri, rows, cols, cnlist, bnlist, iter)

log-sum-product
"""
function log_sum_product!(dec, Lpos, Lpri, dectype, max_iter)
    @unpack rows, cols, cnlist, bnlist = dec

    for iter = 1:max_iter
        if iter == 1
            _init!(Lpri, cols, cnlist)
        else
            _bit_to_check!(Lpos, rows, cols, cnlist)
        end
        _check_to_bit!(dectype, rows, cols, bnlist)
        _calc_posterior_llrs!(Lpos, Lpri, rows, cnlist)
        checksum(Lpos, bnlist)==0 && return 0, iter
    end
    return 1, max_iter
end

# 変数ノード初期メッセージ伝搬
function _init!(Lpri, cols, cnlist)
    mes = 0.0
    for n in eachindex(cnlist) # Threads.@threads
        mes = tanh(Lpri[n]/2.0) # sending message
        for to in cnlist[n]
            cols[n][to] = mes
        end
    end
end


# BitNode -> CheckNode
function _bit_to_check!(Lpos, rows, cols, cnlist)
    # bitnode -> checknode
    Lext = 0.0
    for n in eachindex(cnlist) # Threads.@threads
        for to in cnlist[n]
            Lext = Lpos[n] - rows[to][n]
            cols[n][to] = tanh(Lext/2.0)
            #mes = Lsum - rows[i][n]
            #cols[n][i] = mysign(mes)*phi(mes)
        end
    end
end

# ChekNode -> BitNode
function _check_to_bit!(::LogSumProduct, rows, cols, bnlist)
    mes = 0.0
    for m in eachindex(bnlist) # Threads.@threads
        for to in bnlist[m]
            mes = 1.0
            for from in bnlist[m]
                to == from && continue
                mes *= cols[from][m]
            end
            rows[m][to] = 2.0 * myatanh(mes) # sign(M_j) * min(|M_j|)
        end
    end
end

# ChekNode -> BitNode
function _check_to_bit!(::MinSum, rows, cols, bnlist)
    mes = sign =  val = 0.0
    for m in eachindex(bnlist) # Threads.@threads
        for to in bnlist[m]
            mes  = 1e20
            sign = 1.0
            for from in bnlist[m]
                to == from && continue
                val = cols[from][m]
                sign *= Base.sign(val)
                mes  = ifelse(mes>abs(val), abs(val), mes)
            end
            rows[m][to] = sign * myatanh(mes) # sign(M_j) * min(|M_j|)
        end
    end
end


# Posterior Probablity Calc
function _calc_posterior_llrs!(Lpos, Lpri, rows, cnlist)
    Lsum = 0.0
    for n in eachindex(Lpos)  # Thereads.@threads
        Lsum = Lpri[n]
        for from in cnlist[n] #  チェックノード集合
            Lsum += rows[from][n] # ΣLext + La
        end
        Lpos[n] = Lsum             # Lpri + Lext
    end
end
