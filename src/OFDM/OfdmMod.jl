
mutable struct CP end
mutable struct ZP end
mutable struct Prefix end
mutable struct Postfix end

Vec = Vector # alias
List = Vector
struct IdxsList
    data::List{Vec{Int}}
    pilot::Dict{Int8,List{Vec{Int}}}
end

struct OfdmMod
    n_fft::Int64
    n_gi::Int64
    n_Ts::Int64
    n_dims::Tuple{Int,Int}
    mod_type::ModType
    idxs::IdxsList
    pilot::Union{Nothing,Pilot{<:Layout}}    # 送信パイロット
    gitype::Union{CP,ZP}
    # chest::Dict{String, Symbol}
    dframe::DataFrame
    chest::Dict{Symbol,Union{Symbol,Bool}}
end
Ofdm = OfdmMod

# OfdmMod型のコンストラクタ
function OfdmMod(; n_fft, n_gi, n_dims, mod_type, n_data, pilot=nothing, gitype=:CP, n_null_carriers=0, dc_carrier=false)
    n_Ts = n_fft + n_gi
    n_Tx, n_Rx = n_dims
    # Subcarrier index割当て
    null_idxs = _null_indices(n_null_carriers, n_fft, dc_carrier)
    idxs = setdiff(1:n_fft, null_idxs)
    # idxs = _subcinds(n_fft, n_gi, n_guard_carriers)
    if isnothing(pilot)
        pilot = nothing
    else
        if pilot.df == -1
            pilot.df = div(n_fft, n_gi)
        end
        pilot_idxs = setdiff(1:pilot.df:n_fft, null_idxs)
        pilot = _complete(pilot, pilot_idxs, n_dims[1]) # パイロット型作成
    end
    # Guard Interval
    gitype = gitype == :CP ? CP() : ZP()

    # データフレーム作成
    dframe = DataFrame(eltype(mod_type), n_data, n_fft, n_Tx, pilot, idxs)
    chest = Dict{Symbol,Union{Symbol,Bool}}(
        :flag => isnothing(Pilot),
        :method => :ls,
        :exterp => false,
        :interp => true,
        :average => false,
    )
    OfdmMod(n_fft, n_gi, n_Ts, n_dims, mod_type, idxs, pilot, gitype, dframe, chest)
end


Base.show(io::IO, o::OfdmMod) = print(io, begin
    output = "$(typeof(o))("
    names = fieldnames(typeof(o))
    for name in [:n_fft, :n_gi, :n_dims]
        output *= "$name=>$(getfield(o, name)), "
    end
    output * "...)"
end
)

# 基本パラメタの取得
function getparams(ofdm_t::OfdmMod)
    fnames = (:n_fft, :n_gi, :n_Ts, :n_dims)
    values = getfield.(Ref(ofdm_t), fnames)
    return NamedTuple{fnames}(values)
end

Base.show(io::IO, m::MIME"text/plain", o::OfdmMod) = print(
    io,
    """$(typeof(o)) :
      n_fft    =>  $(o.n_fft)
      n_gi     =>  $(o.n_gi)
      n_dims   =>  $(o.n_dims)
      mod_type =>  $(o.mod_type)
      pilot    =>  $(o.pilot)
  """
)


# nullキャリアインデックスを割り当てる
# ngc分だけ左右均等になるように割り当てる:
# 正の周波数: 0 ~ N/2,
# 負の周波数: N/2+1 ~ N-1]
function _null_indices(ngc, n_fft, dc_carrier=false)
    gidxs = Int[]
    ngc == 0 && return gidxs
    # index割り当て
    i = 0
    while length(gidxs) < ngc
        if (i % 2) == 0
            idx = n_fft - fld(i, 2)
        else
            idx = fld(i, 2) + 1
        end
        push!(gidxs, idx)
        i += 1
    end
    sort!(gidxs)
    return gidxs
end

# パイロットキャリア割り当て
function _pilot_indices(df, guard_idxs, align=:left, style=:equispaced) where {T<:Layout}
    pilot_idxs = Int[]

    if style == :equispaced
        # index計算
        if align == :left
            i = 1
            while 1 <= i <= n_fft
                if i in guard_idxs
                    i += 1
                    continue
                end
                push!(pinds, i)
                i += df
            end
        elseif align == :right
            i = n_fft
            while 1 <= i <= n_fft
                if i in guard_idxs
                    i -= 1
                    continue
                end
                push!(pilot_idxs, i)
                i -= df
            end
        else
            @error ""
        end
    else
        style == :random
        cands = setdiff(1:n_fft, guard_idxs)
        append!(pilot_idxs, sample(cands, pilot.numofpilots, replace=false, ordered=true))
        return pinds
    end

    return sort!(pinds)
end

# OFDMフレームの生成
function gen_frame(ofdm::Ofdm, n_data_sym=0; n_ofdm_sym=0, prealloc=false)
    if n_data_sym == 0 && n_ofdm_sym == 0
        error("datalenもしくはframelenを1以上の値で入力してください.")
    end
    if n_data_sym > 0
        n_ofdm_sym = calc_frame_size(ofdm, ofdm.pilot, n_data_sym)
    end
    n_fft, ntx = ofdm.n_fft, ofdm.n_dims[1]
    dframe = zeros(ComplexF64, n_fft, n_ofdm_sym, ntx)
    idxs = zeros(Int8, n_fft, n_ofdm_sym)
    numofpilot = prealloc ? _pilot_alloc!(ofdm, frame) : 0
    ofdm.dframe.arr = dframe
    ofdm.dframe.idxs = idxs
    ofdm.dframe.on = Dict{Symbol,Int}(:data => 0, :pilot => numofpilot, :null => 0)
end


# フレーム長の計算
function calc_frame_size(ofdm::Ofdm, pilot::Nothing, ndata::Int)
    data_size = div(ndata, ofdm.n_dims[1])
    dinds = ofdm.idxs[:data]
    rem = data_size # 送信データベクトル数
    numofindex = ofdm.n_fft - length(ofdm.idxs[:guard])
    frame_size = 0
    n_data_carrier = length(dinds)
    frame_size = ceil(Int, data_size / n_data_carrier)
end

function calc_frame_size(ofdm::Ofdm, pilot::Pilot, ndata::Int)
    data_size = div(ndata, ofdm.n_dims[1])
    dinds = ofdm.idxs[:data]
    n_state = size(pilot.order, 2)
    numofpilots = length.(pilot.gidxs)
    state = t = 1
    t0, dt = pilot.t0, pilot.dt
    while rem > 0
        numofpilot = 0
        if isalloc(ofdm.pilot, t)
            for n_tx in axes(pilot.order, 1)
                id = pilot.order[n_tx] # get pilot-group id
                id == 0 && continue
                numofpilot += numofpilots[id]
            end
            state = state % n_state + 1
        end
        rem -= numofindex - numofpilot
        t += 1
        frame_size += 1
    end
    return frame_size
end

function calc_frame_size(ofdm::Ofdm, pilot::Pilot{LTE}, ndata::Int)
    data_size = div(ndata, ofdm.n_dims[1])
    didxs = ofdm.idxs[:data]
    ndacar = length(didxs) # 有効サブキャリア数
    nstate = size(pilot.order, 2) # 状態数

    @unpack t0, dt, offset, order = pilot
    j = 1
    prev_j = 1
    state = 1
    frame_size = 0
    rem = data_size
    while rem > 0
        if ((j - t0) % dt) == 0 || (j == prev_j + offset)
            npicar = 0 # パイロットキャリア数
            for k in axes(order, 1)
                id = order[k, state]
                npicar += length(pilot.gidxs[id])
            end
            rem -= ndacar - npicar
            state = (state % nstate) + 1
        else
            rem -= ndacar
        end
        if j == (prev_j + offset)
            prev_j += dt
        end
        j += 1
        frame_size += 1
    end
    return frame_size
end

function getframesize(n_data, n_Tx, pilot::Pilot, subc_idxs)
    # フレーム長計算
    data_size = div(n_data, n_Tx) # 1アンテナストリームのデータシンボル数
    @unpack offset, group_offset, dt, offset, order = pilot

    i = 1
    prev_i = 1
    state = 1
    n_state = size(order, 2)
    frame_len = 0
    n_data_carriers = length(subc_idxs) # データキャリア数
    rem = data_size
    while rem > 0
        n_pilot_carriers = 0 # パイロットキャリア数
        if ((i - offset) % dt) == 0 || (i == prev_i + group_offset)
            # パイロットあり
            for k in axes(order, 1)
                group_idx = order[k, state] # パイロット群のindex
                n_pilot_carriers += length(pilot.group_idxs[group_idx]) # 挿入パイロット数
            end
            state = (state % n_state) + 1
        end
        rem -= n_data_carriers - n_pilot_carriers
        if i == (prev_i + offset)
            prev_i += dt
        end
        i += 1
        frame_len += 1
    end
    return frame_len
end

function getframesize(n_data::Int, n_Tx::Int, pilot::Nothing, subc_idxs)
    n_data_carriers = length(subc_idxs)
    total = frame_len = 0
    while total < n_data
        total += n_data_carriers
        frame_len += 1
    end
    return frame_len
end

# 事前に割り当てるインデックスを決定する
function _indexing!(ofdm::Ofdm, n_data, n_sym)
    n_fft, n_dims, pilot = ofdm.n_fft, ofdm.n_dims, ofdm.pilot
    n_dataset = div(n_data, n_dims[1])
    fidxs = zeros(Int8, n_fft, n_sym)
    rem = n_dataset

    # パイロットなし
    if isnothing(pilot)
        n = 1
        st = 1
        for j in axes(fidxs, 2)
            for i in ofdm.dinds
                rem == 0 && break
                fidxs[i, j] = -1
                rem -= 1
            end
            rem == 0 && break
        end
    else # パイロットあり
        dinds = ofdm.dinds
        order = pilot.order
        n_state = size(order, 2)
        time_table = make_time_table(pilot, size(fidxs, 2))
        state_table = make_state_table(pilot, size(fidxs, 2)) # stateテーブル

        dt = pilot.dt
        t0 = 1 + pilot.offset
        for j in axes(fidxs, 2)
            state = state_table[j]
            id = time_table[j]
            if state > 0
                for k in axes(time_table, 1)
                    id = time_table[k, j]
                    if id > 0
                        for i in pilot.gidxs[id]
                            fidxs[i, j] = k
                        end
                    end
                end
            end
            for i in dinds
                rem == 0 && break
                fidxs[i, j] > 0 && continue
                fidxs[i, j] = -1
                rem -= 1
            end
            rem == 0 && break
        end
    end
end


# パイロットデータ割り当て
function _pilot_alloc!(ofdm::Ofdm, frame)
    isnothing(ofdm.pilot) && return

    # パイロットシンボル生成
    n_pilot = length(ofdm.pinds)
    pilot = ofdm.pilot
    pilot_seqs = gen_random_pilot(pilot)

    # 割り当て
    group_id = 0
    count = 0
    time_table = make_time_table(pilot, size(frame, 2))
    for j in axes(frame, 2) # nsym
        for k in axes(frame, 3) # ntx
            group_id = time_table[k, j] # group id
            if group_id > 0
                for (i, idx) in enumerate(pilot.gidxs[group_id])
                    frame[idx, j, k] = pilot_seqs[group_id][i]
                    frame.idxs[idx, j] = k
                    count += 1
                end
            end
        end
    end
    frame.on[:pilot] = count
    count
end
