using Memoize
using SparseArrays
using UnPack
using DataStructures
ODict = OrderedDict

# BCJR復号用パラメータ
struct BCJRParams
    α::Vector{Vector{Float64}}
    β::Vector{Vector{Float64}}
    γ::Vector{SparseMatrixCSC{Float64,Int64}}
    edges::Dict{Symbol, Dict{Int,Vector{Int}}}  # エッジ
    paths::Dict{Symbol, ODict{Int,ODict{Int,Vector{Tuple{Int,Int}}}}}  # 入力/出力に対応するパスリスト
    trans::OrderedDict{NTuple{2,Int}, NTuple{2,Int}} # 遷移に対応する入出力シンボル
end
function BCJRParams(trellis, N)
    numstates = trellis.nstate
    α = [zeros(numstates) for _ in 1:N+1]
    β = [zeros(numstates) for _ in 1:N+1]
    γ = [spzeros(numstates,numstates) for _ in 1:N]
    edges, paths = _make_edges_and_paths(trellis)
    trans = make_transmat(trellis, dir=:forward)
    BCJRParams(α, β, γ, edges, paths, trans)
end

Base.show(io::IO, m::MIME"text/plain", bcjr::BCJRParams) = print(io,
"""$(typeof(bcjr))$(fieldnames(typeof(bcjr)))""")

# 状態遷移のパス(エッジ)のリスト生成
function _make_edges_and_paths(trellis)
    @unpack nextstates, outputs = trellis
    numstates = trellis.nstate
    numinputs = 2^trellis.size[1]
    numoutputs = 2^trellis.size[2]

    forward  = Dict{Int,Vector{Int}}(k=>Int[] for k in 0:numstates-1)
    backward = Dict{Int,Vector{Int}}(k=>Int[] for k in 0:numstates-1)
    base = ODict{Int,Vector{NTuple{2,Int}}}(k=>Tuple{2,Int}[] for k in 0:1)
    input  = ODict{Int,ODict{Int,Vector{NTuple{2,Int}}}}(ri=>deepcopy(base) for ri in 1:trellis.size[1])
    output = ODict{Int,ODict{Int,Vector{NTuple{2,Int}}}}(ro=>deepcopy(base) for ro in 1:trellis.size[2])
    for n in eachindex(nextstates)
        out = outputs[n] # 出力シンボル
        inp ,ps = divrem(n-1,numstates) # 入力シンボル, 前状態
        ns = nextstates[n]
        push!(forward[ns], ps)
        # inputを保存
        for i in 1:trellis.size[1]
            digit = 1 << (trellis.size[1]-i) # bit桁
            j = (inp&digit) > 0
            push!(input[i][j], (ps,ns))
        end
        # outputを保存
        for i in 1:trellis.size[2]
            digit = 1 << (trellis.size[2]-i)
            j = (out&digit) > 0
            push!(output[i][j], (ps,ns))  # 出力シンボルに対応するペア
        end
    end
    # 後ろ向きの状態推移するノードを状態別に分類
    for i in axes(nextstates,1)
        ps = i - 1 # 前状態
        for j in axes(nextstates,2)
            ns = nextstates[i,j]
            # psとエッジを持つ次状態を保存(入力に依存しない)
            push!(backward[ps], ns)
        end
    end
    edges = Dict(:forward=>forward, :backward=>backward)
    paths = Dict(:input=>input, :output=>output)
    return edges, paths
end

# 入力/出力に対応する遷移パスを辞書として生成する
function _make_paths(nextstates, outputs)
    for i in axes(nextstates,1)
        for j in axes(nextstates,2)
            ns = nextstates[i,j] # 次状態
            input_mat[i,ns+1]  = j
            output_mat[i,ns+1] = outputs[i,j]+1
        end
    end
end

# path: Dict{Int, Vector{Tuple{Int,Int}}}
# edge: NamedTuple{(:forward,:backward), Tuple{Dict{NTuple{N,Int}},Dict{NTuple{N,Int}}} }
#

"""
最大事後確率復号器(APP)オブジェクト
"""
struct APPDecoder <: AbstractDecoder
    trellis::Trellis     # 符号化トレリス
    params::BCJRParams   # BCJRパラメータ
end

abstract type DecType end

struct Viterbi
end

struct APP
end

struct Decoder{T}
    # trellis:Trellis
end


# コンストラクタ
function APPDecoder(trellis::Trellis; numbits, method="max", isterm=true)
    trellis_length =  numbits + trellis.numtailbits # トレリス長
    operator = getoperator(method) # 復号器の作用素
    params = BCJRParams(trellis, trellis_length)
    APPDecoder(trellis, params)
end

Base.show(io::IO, appdec::APPDecoder) = print(io,
"""$(typeof(appdec))""")
Base.show(io::IO, m::MIME"text/plain", appdec::APPDecoder) = print(io,
"""$(typeof(appdec))""")


# 作用素を選択する
function getoperator(method)
    if method=="max"
        return Base.max
    elseif method=="max*"
        error()
    elseif method=="MAP"
        error()
    end
end

abstract type Operator end

struct maxlogMAP <: Operator
end
struct Jacobian <: Operator
end
struct TruleMAP <: Operator
end


"""
受信畳み込み符号に対する最大事後確率復号器
"""
function (dec::APPDecoder)(Lch; target=:input) where T<:AbstractArray
    # 前処理
    trellis = dec.trellis
    Lch = reshape(Lch, trellis.size[2], :)
    trellis_len = size(Lch,2) # トレリスの長さ
    if target == :input
        calc = :input
        bit_size = trellis.size[1]
    else
        calc = :output
        bit_size = trellis.size[2]
    end
    Lpos = Array{Float64}(undef, bit_size, trellis_len) # 出力配列(事後尤度比)

    # BCJR復号の実行
    bcjr = dec.params #
    terminate = isterm(trellis)
    Lpos = decode!(bcjr, Lpos, Lch, calc, terminate) # BCJRアルゴリズムによるMAP復号
    if target == :input
        return @view Lpos[1:end-trellis.numtailbits]
    else
        return vec(Lpos)
    end
end



"""
 BCJRアルゴリズムによる最大事後確率復号
"""
function decode!(params::BCJRParams, Lpos, Lch, calc, terminate; Lpri=nothing)
    @unpack α, β, γ, edges, paths, trans = params
    _init!(Lch, Lpri, α, β, γ, edges, paths, terminate)
    _forward!(Lch, Lpri, α, γ, edges, paths, trans, max)
    _backward!(Lpos, α, γ, β, edges, paths, calc, max)
    return Lpos
end

# α, β　初期化
function _init!(Lch, Lpri, α, β, γ, edges, paths, terminate)
    α[1][1] = 0.0; α[1][2:end] .= -Inf
    if terminate
        β[end][1] = 0.0; β[end][2:end] .= -Inf
    else
        β[end][:] .= -log(length(β[end]))
    end
end

# γ(チャネル尤度の計算)
function _calc_gamma(Lch, k)
    lp = 0.0 # 対数尤度
    for i in reverse(eachindex(Lch))
        digit = 1 << (i-1)
        x = (k&digit) > 0 # binary
        lp += Lch[end-i+1]/2.0 * (2x-1)
    end
    lp
end


# 前向き計算
function _forward!(Lch, Lpri, α, γ, edges, paths, trans, operator)
    N = size(Lch,2) # トレリス長
    numstates = length(α[1]) # 状態数
    list = edges[:forward] # 前向き遷移のパス

    for n = 1:N
        # γの更新
        for ((i,j), inout) in trans
            @views γ[n][i+1,j+1] = _calc_gamma(Lch[:,n], inout[2])
        end
        if !isnothing(Lpri) && (size(Lpri,2) <= n)
            for ((i,j), inout) in trans
                @views γ[n][i+1,j+1] += _calc_gamma(Lch[:,n], inout[1])
            end
        end

        # αの更新
        αmax = -Inf
        for (j,i_idxs) in list
            tmp = -Inf
            for i in i_idxs
                tmp = operator(tmp, γ[n][i+1,j+1] + α[n][i+1])
            end
            α[n+1][j+1] = tmp
            αmax = ifelse(tmp>αmax, tmp, αmax) # αの最大値を保存
        end
        for j in 0:numstates-1
            α[n+1][j+1] -= αmax
        end
    end
end

# 後ろ向き計算
function _backward!(Lapp, α, γ, β, edges, paths, calc, operator)
    N = size(Lapp,2)
    numstates = length(α[1])
    list = edges[:backward]
    path = paths[calc]

    for n in N:-1:1
        # 事後尤度比の計算
        @views _calcapp!(Lapp[:,n], α[n], γ[n], β[n+1], path, operator)
        n == 1 && continue

        βmax = -Inf
        for (i,idxs) in list
            tmp = -Inf
            for j in idxs
                tmp = operator(tmp, γ[n][i+1,j+1]+β[n+1][j+1])
            end
            β[n][i+1] = tmp
            βmax  = ifelse(tmp>βmax, tmp, βmax) # αの最大値を保存
        end
        for i in 0:numstates-1
            β[n][i+1] -= βmax
        end
    end
end

# k番目bitに対する事後LLRを計算
function _calcapp!(lpos, α, γ, β, path, operator)
    # Posterior LLR at k-th bit
    for k in eachindex(posterior_llr) 
        den = num = -Inf
        for (i,j) in path[k][1]
            # P(x = +1)
            den = operator(A, α[i+1]+γ[i+1,j+1]+β[j+1])
        end
        for (i,j) in path[k][0]
            # P(x = -1)
            num = operator(B, α[i+1]+γ[i+1,j+1]+β[j+1])
        end
        lpos[k] = den - num # Lapp ~ max{p(x=+1|y)} - max{p(x=-1|y)}
    end
end
