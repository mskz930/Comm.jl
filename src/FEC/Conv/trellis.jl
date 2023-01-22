
"""
トレリス多項式
"""
struct Polynomial{T} <: AbstractVector{T}
    v::Vector{T}
    Polynomial(v::Vector{T}) where T= Polynomial{T}(v)
end

Base.length(p::Polynomial) = length(p.v)
Base.size(p::Polynomial) = (length(p.v),)
Base.getindex(p::Polynomial, i...) = p.v[i...]


function Base.convert(::Type{Polynomial}, vs::Vector{String}; to=2, from=8)
    vs = parse.(Int, vs, base=from) # n進数を10進数に変換
    vs = digits.(vs, base=to) # 2進数
    vs = reverse.(vs)
    [Polynomial(v) for v in vs]
end


"""
トレリス

Example: 
    poly2trellis()
"""
struct Trellis
	n_inputs::Int
	n_outputs::Int
    size::Tuple{Int64,Int64}    # 入出力次元
    nstates::Int64              # 状態数
    ntap::Int64                 # タップ数
    K::Vector{Int64}            # 拘束長
    numtailbits::Int64          # 終端ビット長
    nextstates::Matrix{Int64}   # 次状態リスト
    outputs::Matrix{Int64}      # 出力リスト
    terminate::Bool             # トレリス終端
end

show(io::IO, t::Trellis) = print(io, """$(t.nstates)-states Trellis""")
show(io::IO, m::MIME"text/plain", t::Trellis) = print(io, """$(t.nstates)-states Trellis""")

# polynomial(8進数)からトレリスオブジェクトを生成
function poly2trellis(poly::AbstractArray; fbpoly=Int[], terminate=true, puncture=false)
	ndims(poly) == 1 && (poly = reshape(poly, 1, length(poly)))
	# fbpoly: feedback polynomial
    input_size, output_size = size(poly)
    if poly isa AbstractVector
        poly = poly'
    end

    poly, fbpoly = string.(poly), string.(fbpoly)
    poly = parse.(Int, poly, base=8)
    fbpoly = parse.(Int, fbpoly, base=8)
    mask = isempty(fbpoly) ? poly .< 0 : poly .== fbpoly 
    K = maximum(ndigits.(poly, base=2), dims=2)[:] # 拘束長(constraint length)
    numtailbits = terminate ? (maximum(K)-1) : 0 # 終端ビット数
    ntap = sum(K) - length(K) # 状態数
    nstates = 2^ntap

    # 状態遷移行列の取得
    isfedback = !isempty(fbpoly)
    nextstates, outputs = _make_nextstates(poly, fbpoly, mask, K, isfedback)
    return Trellis(
            input_size, 
            output_size, 
            (input_size, output_size), 
            nstates, 
            ntap, 
            K, 
            numtailbits, 
            nextstates, 
            outputs, 
            terminate)
end

rate(t::Trellis) = t.n_inputs // t.n_outputs
is_terminated(t::Trellis) = t.terminate # トレリス終端判定


# 出力符号語ビット数を返す
function num_of_outputs(trellis, nbits)
    ceil(Int, nbits/float(getrate(trellis))) + trellis.numtailbits*trellis.size[2]
end

# 前状態の保存
function _make_prevstates(nextstates)
    prevstates = fill(-1, size(nextstates))
    for i in axes(prevstates,1)
        ps = i - 1 # 前状態
        for j in axes(prevstates,2)
            ns = nextstates[i,j] # 次状態
            prevstates[ns+1,j] = ps # 前状態
        end
    end
    prevstates
end


# 次状態配列を作成する
function _make_nextstates(poly, fbpoly, polymask, K, isfeedback)
    # パラメータ
    insize, outsize = size(poly)
    L = sum(K) - length(K) # 総レジスタ数
    n_state = 2^L; n_input_states = 2^insize # 状態数

    # 出力配列確保
    nextstates = zeros(Int64, 2^L, 2^insize) # 次状態配列
    outputs = zeros(Int64, 2^L, 2^insize) # 出力配列

    # 現状態(i)についての次状態、出力を求める。
    for m in 0:n_state-1 #
        fbval = 0 # フィードバック値
        if isfeedback
            for i in insize
                temp = m >> (L - K[end-i+1] + 1)
                fbval += _xor_sum(fbpoly[i], temp) << (insize-i)
            end
        end
        for n in 0:n_input_states-1
            next_state = output = 0
            inp = dec2bin(n ⊻ fbval; pad=insize)
            # 入力シンボルに対する状態遷移
            for i in 1:insize
                temp = m >> (L - K[end] + 1)  # 次状態
                register = temp + (inp[i] << (K[i]-1)) # レジスタ状態
                for j in 1:outsize # 出力計算
                    output += _xor_sum(register, poly[i,j]) << (outsize-j) # XOR_SUM
                end # for
                register >>= 1 # register shift
                next_state += register << (L - sum(K[end-i+1:end]) + i)
            end # for
            nextstates[m+1, n+1] = next_state
            outputs[m+1, n+1] = output
        end # for

    end # for
    nextstates, outputs
end # function

# XOR_SUM
function _xor_sum(x, y)
    z = x & y # XOR
    s = 0 # xor_sum
    while z>0
        s ⊻= z & 1 # AND
        z >>= 1 # bit shift
    end
    s
end

# 配列を任意のサブ配列(spの要素おきに)に分割
function _split(arr::AbstractVector, sp::Union{Integer, AbstractArray}; dir="left")
    @assert sum(sp)==length(arr) # 分割の和が配列要素を超える場合
    if isa(sp, Integer)
        return [arr]
    else
        subarr = Array{Array{eltype(arr),1},1}(undef,length(sp)) # 出力配列
        if dir == "left"
            pos = 1
            for i in 1:length(sp)
                subarr[i] = arr[pos:pos+sp[i]-1]
                pos += sp[i]
            end
        else
            pos = length(arr)
            for i in 1:length(sp)
                subarr[i] = arr[pos-sp[i]+1:pos]
                pos -= sp[i]
            end
        end
        return subarr
    end
end

# 状態遷移行列の生成
function make_transmat(trellis::Trellis; dir=:forward)
    @unpack nstates, outputs, nextstates = trellis
    transmat = OrderedDict{Tuple{Int,Int}, Tuple{Int,Int}}()

    m,n = size(nextstates)
    if dir == :forward
        for ps in 0:m-1
            for inp in 0:n-1
                ns = nextstates[ps+1,inp+1]
                out = outputs[ps+1,inp+1]
                transmat[(ps,ns)] = (inp,out)
            end
        end
    else
        for ps in 0:m-1
            for inp in 0:n-1
                ns = nextstates[ps+1,inp+1]
                out = outputs[ps+1,inp+1]
                transmat[(ps,ns)] = (inp,out)
            end
        end
    end
    sort!(transmat, by=x->x[1])
end
