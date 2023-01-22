
# Ofdm Frame Type
mutable struct Frame
    array::Array{ComplexF64,3}
    subcind::Array{Int64,2}
    pilot_is_inserted::Bool
    ndata::Vector{Int64}
end


function Frame(ofdm, ndata=0; dataseq=[], framelen=0, mapping=true, pilot_insert=false)
    if ndata > 0 && framelen == 0
        framelen = calc_framesize(ofdm, ndata)
    end
    frame_array  = zeros(ComplexF64, ofdm.nfft, framelen, ofdm.ndim[1])
    subcind      = zeros(Int, ofdm.nfft, framelen)

    mapping      && _make_subcind!(ofdm, subcind, ndata)
    pilot_insert && _pilot_alloc!(ofdm, array, subcind)

    Frame(frame_array, subcind, pilot_insert, [0])
end

function crete_frame(ofdm::Ofdm, n_data, framelen=0)
    if ndata > 0 && framelen == 0
        framelen = calc_framesize(ofdm, ndata)
    end
    frame  = zeros(ComplexF64, ofdm.nfft, framelen, ofdm.ndim[1])
    if ofdm.frame.pilot_is_preallocate
        _pilot_alloc!()
    end
    return frame
end


# メソッド
Base.length(frame::Frame)  = length(frame.array)
Base.size(frame::Frame)    = size(frame.array)
Base.size(frame::Frame, x) = size(frame.array, x)
Base.axes(frame::Frame, x) = axes(frame.array, x)
Base.getindex(frame::Frame, inds...) = getindex(frame.array, inds...)
Base.setindex!(frame::Frame, v, inds...) = setindex!(frame.array, v, inds...)
Base.copy(frame::Frame) = Frame((getfield(frame, k) for k in fieldnames(Frame))...)

# フレーム長の計算
function calc_framesize(ofdm::Ofdm, n_data)
    n_dataset = div(n_data, ofdm.ndim[1])
    pilot = ofdm.pilot
    rem = n_dataset # 送信データセット数
    framesize = 0

    # パイロット挿入なしの場合
    if isnothing(pilot)
        n_alloc_data = length(ofdm.dinds)
        while rem > 0
            framesize += 1
            rem -= n_alloc_data
        end
    else
        # パイロット挿入する場合
        n_state = size(pilot.order,2)         # パターン数
        numofpilots  = length.(pilot.inds)    # 各送信パターンのパイロットシンボル数
        max_data_len = length(ofdm.dinds)     # 最大データシンボル数
        state = t = 1
        while rem > 0
            numofdata = max_data_len
            framesize += 1
            if (t - pilot.t0) % pilot.dt == 0
                for n_tx in axes(pilot.order,1)
                    id = pilot.order[n_tx] # pilot group id
                    id == 0 && continue
                    numofdata -= numofpilots[id]
                end
                state = state % n_state + 1
            end
            rem -= numofdata
            t += 1
        end
    end
    framesize
end

# subcindに割り当てindexを書き込む
function _make_subcind!(ofdm::Ofdm, subcind, n_data)
    nfft, ndim, pilot   = ofdm.nfft, ofdm.ndim, ofdm.pilot
    n_dataset = div(n_data, ndim[1])
    rem = n_dataset

    if isnothing(pilot) # パイロットなし
        n = 1; st= 1
        for j in axes(subcind,2)
            for i in ofdm.dinds
                rem ==0 && break
                subcind[i,j] = -1
                rem -= 1
            end
            rem == 0 && break
        end
    else # パイロットあり
        order  = pilot.order
        n_state      = size(order,2)
        time_table   = make_time_table(pilot, size(subcind,2))
        state_table  = make_state_table(pilot, size(subcind,2)) # stateテーブル
        dinds = [setdiff(ofdm.inds, s) for s in pilot.inds]
        t0 = pilot.t0
        dt = pilot.dt
        for j in axes(subcind,2)
            state = state_table[j]
            id = time_table[j]
            if state > 0
                for k in axes(time_table,1)
                    id = time_table[k,j]
                    if id > 0
                        pinds = pilot.inds[id]
                        for i in pinds
                            subcind[i,j] = k # tx antenna
                        end
                    end
                end
            end
            # データシンボルindexの配置
            for i in ofdm.dinds
                rem==0 && break
                subcind[i,j] = -1
                rem -= 1
            end
            rem==0 && break
        end
    end
end


# パイロットデータ割り当て
function _pilot_alloc!(ofdm::Ofdm, frame::Frame)
    pilot = ofdm.pilot
    pilot==nothing && return

    n_pilot = length(ofdm.inds[:pilot])
    ofdm.subcind[:pilot] = Vector{Int}[]

    # パイロットシンボル生成
    time_table = _make_time_table(pilot, size(frame,3), size(frame,2))
    id = 0
    for j in axes(frame,2)
        tmp = Int[]
        for k in axes(frame,3)
            id = time_table[j,k] # group id
            if id > 0
                for (idx,i) in enumerate(pilot.inds[id])
                    frame[i,j,k] = pilot_symbols[id][idx]
                end
                append!(tmp, pilot.inds[id])
            end
        end
        push!(ofdm.subcind[:pilot], sort!(tmp))
    end
    frame.is_pilot_inserted = true
end

# パイロットマッピング
function pilot_map!(ofdm::Ofdm, frame, subcind)
    n_pilot = length(ofdm.inds[:pilot])
    n_pilot==0 && return
    pilot = ofdm.pilot
    pilot_symbols = rand(Qam(4,isnorm=false), n_pilot)
    pilot_seq = _split(pilot_symbols, length(pilot.inds))
    pilot.symbols = pilot_symbols

    # パイロットマッピング
    (size(pilot.order, 1) < size(frame, 3)) && error("pilot.orderに誤りがあります.")
    n_state = size(pilot.order, 2)
    state = 1
    j = pilot.t0 # 初期挿入index
    while j <= size(frame, 2) # 時間軸
        # パイロット挿入
        for k in axes(frame, 3) # 空間軸
            id = pilot.order[k,state] # pilot group ID
            id == 0 && continue
            for (l, i) in enumerate(pilot.inds[id])
                frame[i,j,k] = pilot_seq[id][l]
            end
        end
        j += pilot.dt # 次点の挿入index
        state = state % n_state + 1 # 次状態
    end
    frame
end


# パイロット挿入時間id列の取得
function _get_pilot_time_inds(pilot::Pilot, frame_size)
    n_tx = length(pilot.order)
    time_inds = zeros(Int64, frame_size, n_tx) # time index
    for j in 1:n_tx
        p = 1;
        for i in pilot.init_idx:pilot.dt:frame_size
            time_inds[i,j] = pilot.order[j][p]
            p = p%pilot.n_groups + 1
        end
    end
    time_inds
end

# index
function get_inds(frame::Frame, target=:data, id=0)
    subcind = frame.subcind
    if target==:data
        return findall(x->x<0, subcind)
    elseif target==:pilot
        if id > 0
            findall(x->x==id, subcind)
        elseif id == 0
            findall(x->x>0, subcind[:,j])
        end
    elseif target==:data_pilot
        findall(!iszero, subcind[:,j])
    else
        findall(iszero, subcind[:,j])
    end
end

# subcindからCartesianIndex列を作成する
function get_inds_list(frame::Frame, target=:data, id=0)
    subcind = frame.subcind
    if target==:data
        return [findall(x->x<0, subcind[:,j]) for j in axes(subcind,2)]
    elseif target==:pilot
        if id > 0
            [findall(x->x==id, subcind[:,j]) for j in axes(subcind,2)]
        elseif id == 0
            return [findall(x->x>0, subcind[:,j]) for j in axes(subcind,2)]
        end
    else
        return [findall(x->x==0, subcind[:,j]) for j in axes(subcind,2)]
    end
end
