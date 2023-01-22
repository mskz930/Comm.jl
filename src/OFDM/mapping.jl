# mapping.jl

function alloc_indices(ofdm::OfdmMod)
    dframe = ofdm
    frame = dframe.A
    tab = dframe.table_index
    pilot = ofdm.pilot
    
    for j = axes(frame,2)
        pinds = []
        # パイロット割り当て
        gid = gids[j]
        if gid > 0
            pinds = pilot.indices[gid]
            for idx in pinds
                tab[i,j] = gid
            end
        end

        # データ割り当て
        for i in ofdm.idxs.data
            if !(i in pinds)
                tab[i,j] = -1
            end
        end
    end
    
end

# サブキャリア割当て
function mapping(ofdm::Ofdm, tx_data; allocation=:static)
    tx_frame = copy(ofdm.dframe.A)
    # _init_frame_idxs!(ofdm.dframe)
    tx_frame = mapping!(tx_frame, tx_data, ofdm.dframe, ofdm.pilot)
    # dataindex, pilot indexを割り当て直す
    # allocation == :dynamic && _pilot_realloc!(ofdm.pilot, ofdm.idxs[:data], ofdm.idxs[:pilot])
    # ofdm.pilot.symbols = pilot_seq = gen_random_pilot(ofdm.pilot)
    # pilot_seqs = _split(pilot_seq, ofdm.pilot.n_group)
    # tx_frame = mapping!(ofdm, tx_frame, data_seq, pilot_seqs)
end

function mapping!(tx_frame, dframe::DataFrame, pilot::Pilot)
    
    for j = axes(tx_frame,2)
        for i = axes(tx_frame,1)
        end
    end
end


# サブキャリアにデータシンボルを割り当てる
function mapping!(tx_frame, tx_data, dframe::DataFrame, pilot::Nothing)
    n_data  = length(tx_data)
    tx_data = serial_to_parallel(tx_data, size(tx_frame, 3), :elm)
    data_idxs_list = dframe.idxs_list[:data]
    idx_table = dframe.idx_table

    n_dummy = 0
    count = 0
    for j in axes(tx_frame,2) # time
        push!(data_idxs_list, Int8[])
        for i in data_idxs # freq
            if count <= size(tx_data, 2)
                for k in axes(tx_frame,3) # space
                    tx_frame[i,j,k] = tx_data[k,n]
                end
                push!(data_idxs_list[j], i)
                idx_table[i,j] = -1
                count += 1
            else
                # ダミーデータの挿入
                for k in axes(tx_frame,3) # space
                    tx_frame[i,j,k] = tx_data[(i+k) % length(tx_data)+1]
                end
                idx_table[i,j] = -2
                n_dummy += 1
            end
        end
        count >= size(tx_data,2) && break
    end
    n_data, n_dummy
end


# サブキャリアにデータシンボルおよびパイロットシンボルを割り当てる
function mapping!(tx_frame, tx_data, dframe::DataFrame, pilot::Pilot)
    ndata = length(tx_data)
    ndummy = ndata % size(tx_frame,3)
    data_idxs = ofdm.idxs[:data] # 有効なデータindex

    tx_data = serial_to_parallel(tx_data, size(tx_frame,3), :random)
    pcount = dcount = 0
    dind = j = st = 1
    time_table = make_time_table(pilot, size(tx_frame,2))
    tidxs = sum(time_table, dims=1) # pilot time slot
    for j in axes(tx_frame,2) # time
        # パイロット割り当て
        if tidxs[j] > 0
            for k in axes(tx_frame,3) # space
                id = time_table[k,j]
                id == 0 && continue
                for (i, val) in zip(pilot.gidxs[id], pilot_seqs[id])
                    tx_frame[i,j,k] = val
                    ofdm.dframe.idxs[i,j] = k
                    pcount += 1
                end
            end
        end
        # データ割り当て
        for i in data_idxs
            ofdm.dframe.idxs[i,j] > 0 && continue
            if dind <= size(tx_data,2)
                # data symbol mapping
                for k in axes(tx_frame,3) # space
                    tx_frame[i,j,k] = tx_data[k,dind]
                end
                ofdm.dframe.idxs[i,j] = -1
            else
                # dummy symbol mapping
                for k in axes(tx_frame,3) # space
                    tx_frame[i,j,k] = tx_data[(i * k) % length(tx_data) + 1]
                    ndummy += 1
                end
                ofdm.dframe.idxs[i,j] = -2
            end

            dind += 1 # next dataset
        end
        dind > size(tx_data,2)
    end # while
    dcount = dind*size(tx_frame,3)
    ofdm.dframe.on[:data] = ndata
    ofdm.dframe.on[:dummy] = ndummy
    return tx_frame
end


# データシンボルとパイロットシンボルをシャッフルして
# ランダムに割り当てなおす
function _pilot_realloc!(pilot, dinds, pinds)
  nd, np = length(dinds), length(pinds)
  inds = union(dinds, pinds) |> shuffle!
  pinds .= @views sort!(inds[1:np])
  # dinds .= @views sort!(inds[np+1:end])
  pilot.gidxs = _split(pinds, pilot.n_group)
  return pilot
end


"""
サブキャリアデマッピング
"""
function demapping(ofdm::Ofdm, rx_frame)
    ndata = ofdm.dframe.on[:data] # 有効データシンボル数
    idxs = get_idxs(ofdm.dframe, :data)
    rx_data = @views rx_frame[idxs, :]
    return @views rx_data[1:ndata]
end

"""
    data_mapping(params::Ofdm, sig_frame, data)

データシンボルを信号フレームにマッピング
"""
function data_mapping(params::Ofdm, sig_frame::AbstractArray{T,2}, data) where T
    ntx = params.ndim[1]
    data_idxs = params.frame.idxs_list[:data] # subcarrier indices
    cnt = 1

    for j = axes(sig_frame, 2) # time
        for i = data_idxs[j] # freq x space
            for k = (i-1)*ntx+1:i*ntx
                sig_frame[k,j] = data[cnt]
                cnt += 1
            end
        end
    end
    sig_frame
end

function data_mapping(params::Ofdm, sig_frame::AbstractArray{T,3}, data) where T
    data_idxs = params.frame.idxs_list[:data]
    cnt = 1
    for j = axes(sig_frame, 2) # time
        for i = data_idxs[j] # freq
            for k = axes(sig_frame,3) # space
                sig_frame[i,j,k] = data[cnt]
                cnt += 1
            end
        end
    end
    sig_frame
end

# stack信号フレームにデータを挿入する
function data_stacking(params::Ofdm, sig_frame::AbstractArray{T,2}, data) where T
    ntx, n_fft = params.ndim[1], params.n_fft
    data_idxs = params.frame.idxs_list[:data] # subcarrier indices
    cnt = 1

    for j = axes(sig_frame, 2) # time
        for i = data_idxs[j] # freq x space
            for k = i:n_fft:n_fft*ntx
                sig_frame[k,j] = data[cnt]
                cnt += 1
            end
        end
    end
    sig_frame
end