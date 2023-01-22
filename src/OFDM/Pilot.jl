# pilot.jl
abstract type Layout end
abstract type Comb <: Layout end
abstract type Block <: Layout end
abstract type Scatter <: Layout end
abstract type LTE <: Layout end
Layouts = Union{Comb, Block, Scatter, LTE}

"""
Pilot型

example):
    Pilot(layout=:comb, t0=1, dt=2, df=8)

style: :equispaced, :random
"""
mutable struct Pilot{T <: Layout}
    dt::Int                          # 時間周期
    df::Int                          # 周波数間隔
    offset::Union{Nothing,Int}       # グループ内時間オフセット
    group_offset::Int                # グループ内オフセット時間
    isnorm::Bool                     # 電力正規化(true/false)
    order::Matrix{Int}               # 送信順序配列
    idxs::Vector{Int}                # 割り当てられたindex全体
    group_idxs::Vector{Vector{Int}}  # グループごとのindex
    n_group                          # グループ数
    function Pilot{T}(; dt=1, df=-1, offset=0, group_offset=dt, isnorm=false) where T <: Layout
        new{T}(dt, df, offset, group_offset, isnorm)
    end
end
function _complete(pilot::Pilot{T}, pilot_idxs, n_tx) where T
    pilot.idxs = pilot_idxs
    pilot.order = _get_pilot_order(T, n_tx)
    pilot.group_idxs = _to_subgroups(T, pilot.idxs, n_tx)
    pilot.n_group = length(pilot.group_idxs)
    return pilot
end




function Base.show(io::IO, m::MIME"text/plain", pilot::Pilot)
    fields = [:dt, :df, :offset]
    str = "$(typeof(pilot)): \n"
    for field in fields
        str *= "  :$field => $(getfield(pilot, field)) \n"
    end
    print(io, str)
end

# 送信次元を取得
getdim(p::Pilot) = size(p.order,1)

# パイロット配置型を生成
function _get_layout_type(layout::AbstractString)
    layout = Symbol(lowercase(layout))
    if layout == :comb
        return Comb
    elseif layout == :block
        return Block
    elseif layout == :scatter
        return Scatter
    elseif layout == :lte
        return LTE
    else
        error(layout, "は見つかりません。")
    end
end





"""
パイロット送信順序を決定する
   Example)  (Nt,Nr) = (2,2)のとき
    行, 列 : 送信アンテナ次元, 送信周期
    配列値 : 送信パイロットグループ
    OrderMatrix = [1, 2;
                   2, 1]
"""
function _get_pilot_order(layout, n_tx)
    if layout <: Comb
        if n_tx == 1
            order = [1;][:,:]
        elseif n_tx == 2
            order = [1 2;
                     2 1]
        elseif n_tx == 3
            order = [1 0 2 1 0 2;
                     0 2 1 0 2 1;
                     2 1 0 2 1 0]
        elseif n_tx == 4
            order = [1 0 0 2;
                     0 1 0 0;
                     0 2 1 0;
                     0 0 2 1]
        else
            error("n_dims>4はサポートされていません。")
        end

    elseif layout <: Block
        if n_tx == 1
            order=[1;][:,:]
        else
            order=[1 2;
                   2 1]
        end
    elseif layout <: Scatter
        error("準備中です。")

    elseif layout <: LTE
        if n_tx == 1
            order = [1 2]
        elseif n_tx == 2
            order = [1 2;
                     2 1]
        elseif n_tx == 3
            order = [1 0 2 1 0 2;
                     0 2 1 0 2 1;
                     2 1 0 2 1 0]
        elseif n_tx == 4
            order = [1 0 2 0;
                     2 0 1 0;
                     0 1 0 2;
                     0 2 0 1]
        else
            error("n_dims>4はサポートされていません。")
        end
    else
        order =  Matrix{Int}(undef, 0, 0)
    end
    order
end


# layoutに従ってサブグループに分割する
# Vector{Int}で返す
function _to_subgroups(layout, pinds, ntx)
    if (layout <: Comb || layout <: Block) && ntx==1
        n_group = 1
    elseif (layout <: Comb && ntx>1) || layout <: LTE
        n_group= 2
    elseif layout <: Scatter
        n_group = 4
    else
        error(layout, "は存在しません。")
    end
    return _split(pinds, n_group)
end


# ランダムパイロットシンボルの生成
function gen_random_pilot(p::Pilot)
    n_pilot = sum(length(elm) for elm in p.idxs)
    pilot_seq = rand(QamMod(M=4, isnorm=p.isnorm), n_pilot)
    pilot_seq
end


# パイロット周期の計算
function get_period(pilot::Pilot)
    inds =  @views findall(x->x==1, pilot.order[1,:])
    if length(inds) == 1
        period = size(pilot.order,2) * pilot.dt
    else
        period = (inds[2] - inds[1]) * pilot.dt
    end
end

# 挿入timeslotのindexリストを返す
function _get_time_inds(pilot::Pilot, len)
    inds = Int[] # 出力リスト
    i = pilot.t0
    while i <= len
        push!(inds, i)
        i += pilot.dt
    end
    inds
end

# パイロットのtimeslotにおける送信テーブルを作成
function make_time_table(pilot::Pilot, len)
    table   = zeros(Int, size(pilot.order,1), len)
    n_state = size(pilot.order,2)
    j = 1        # state index
    n = pilot.t0 # time-index 初期値
    while n <= len
        table[:,n] .= @view pilot.order[:,j]
        j = j % n_state + 1
        n += pilot.dt + 1
    end
    table
end

# パイロット挿入のタイムテーブル作成
function make_time_table(pilot::Pilot{LTE}, len)
    table   = zeros(Int, size(pilot.order,1), len)
    @unpack t0, dt, offset = pilot
    nstate = size(pilot.order,2)

    state = 1
    for i = t0:dt:len
        table[:,i] .= @view pilot.order[:,state]
        table[:,i+offset] .= @view pilot.order[:,state+1]
        state = mod(state+1, nstate) + 1
    end

    #=
    # pilot.order行列(行: 送信アンテナ, 列: 送信状態)に基づいて挿入timeスロッを決定する
    while n <= len
        if (n+1) <= len
            for i in 1:2
                table[:,n] .= @view pilot.order[:,j]
                j = j % nstate + 1 # 状態(tableの列)
                n += 1
            end
        else
            break
        end
        n += pilot.dt
    end
    =#
    return table
end

# 状態時系列の生成
function make_state_table(pilot::Pilot, len)
    table   = zeros(Int, len)
    n_state = size(pilot.order,2)
    j = 1; k = pilot.t0
    while k <= len
        table[k] = j
        j = j % n_state + 1
        k += pilot.dt
    end
    table
end

# パイロット挿入時間id列の取得
get_time_idxs(p::Pilot, len) = [p.t0:p.dt:len;]

function get_time_idxs(p::Pilot, len::Int)
    @unpack t0, dt = p
    idxs = [t0:dt:len; ]
    return idxs
end
function get_time_idxs(p::Pilot{LTE}, len::Int)
    idxs = zeros(Int, len)
    offset = p.offset
    for state = 1:p.n_group
        group_offset = p.group_offset * (state-1)
        step = p.dt * p.n_group
        for i = 1+offset+group_offset:step:len
            idxs[i] != 0 && error("$(idxs[i])はすでに割り当てられています。")
            idxs[i] = state
        end
    end
    return idxs
end
function get_time_idxs(p::Pilot, len::Int, tx_idx::Int, group_idx::Int)
    idxs = Int[] # time index
    s = 1 # state
    ns = size(p.order,2) # num of states
    for i in p.t0:p.dt:framelen
        p.order[tx_idx, s]==group_idx && push!(idxs, i)
        s = s % ns + 1 # next state
    end
    return idxs
end
function get_time_idxs(p::Pilot{LTE}, len::Int, k::Int, gid::Int)
    idxs = Int[] # time index
    state = 1 # state
    ns = size(p.order,2) # num of states
    for i in p.t0:p.dt:len
        p.order[k, state]==gid && push!(idxs, i)
        state = (state % ns) + 1
        p.order[k, state]==gid && push!(idxs, i+p.offset)
        state = (state % ns) + 1
    end
    return idxs
end

#=
function get_time_idxs(p::Pilot{LTE}, len, tx_idx, group_idx)
    idxs = Int[] # time index
    i = tx_idx
    j = 1
    id = group_idx
    nstate = size(p.order,2) # num of states
    n = p.t0
    while n <= len
        for k in 1:2
            ((n+1)<=len) && id==p.order[i,j] && push!(idxs, n)
            j = j % nstate + 1 # next state
            n += 1
        end
        n += p.dt
    end
    idxs
end
=#

# time indexに対してパイロットの挿入可否を返す
function isalloc(p::Pilot{LTE}, tidx)
    # 先頭パイロットの周期 dt + 1
    tidx < p.t0 && return false
    state = (tidx-p.t0) % (p.dt+2)
    if state == 0 || ((tidx-1) > 0 && state == 1)
        return true
    else
        return false
    end
end

# 最小のコヒーレントパイロット挿入間隔を得る
function mincw(len, n_fft)
    @assert 0 < len <= n_fft
    df = 0 # パイロット挿入間隔
    k = 0; max_l = n_fft

    while true
        k += 1
        max_l = fld(n_fft, k) # 最大周期
        if max_l >= len
            df = k
        else
            break
        end
    end
    return df
end
