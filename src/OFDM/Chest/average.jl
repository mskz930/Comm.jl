
# パイロット平均
function averaging(pilot::Pilot, frame, dir=:time)
    if dir == :time
        averaged = averaging(pilot, frame, Time())
    elseif dir == :freq
        averaged = averaging(pilot, frame, Freq())
    elseif dir == :time_freq

    end
    if size(averaged, 2) == 1
        averaged = dropdims(averaged, dims=2)
    end
    averaged
end

function domain()
    if dir == :time
        Time()
    elseif dir == :freq
        Freq()
    elseif dir == :time_freq
        TimeFreq()
    else
        error()
    end
end

function averaging(pilot::Pilot, rxframe::AbstractArray{T,3}, ::Time) where T
    ntx, nrx = getdim(pilot), size(rxframe,3)
    ngroup = pilot.n_group
    averaged = zeros(T, size(rxframe,1), nrx, ntx)

    # 時間方向平均
    for p in 1:ntx
        for g in 1:ngroup
            tidxs = get_time_idxs(pilot, size(rxframe,2), p, g)
            fidxs = pilot.gidxs[g]
            for q in 1:nrx
                @views mean!(averaged[fidxs,q,p], rxframe[fidxs,tidxs,q])
            end
        end
    end
    averaged
end
function averaging(pilot::Pilot, CFRs::AbstractArray{T,4}, ::Time) where T
    ngroup = pilot.n_group
    averaged = zeros(T, size(CFRs,1), 1, size(CFRs,3), size(CFRs,4))
    # 時間方向平均
    for p in axes(CFRs,4)
        for gid in 1:ngroup
            tidxs = get_time_idxs(pilot, size(CFRs,2), p, gid)
            fidxs = pilot.gidxs[gid]
            for q in axes(CFRs,3)
                @views mean!(averaged[fidxs,1,q,p], CFRs[fidxs,tidxs,q,p])
            end
        end
    end
    averaged
end


# 平均化
function _average(v::AbstractArray, ref)
    sum(view(v,ref)) / length(ref)
end


# 平均化関数
function _average(v::AbstractArray{T}, ind::AbstractVector, w_size::Integer) where T
    p = copy(v)
    wh = Int(w_size-1/2) # half size
    for j in axes(v,2)
        for i in axes(v,1)
            temp = zero(T)
            horz=ramp(j-wh):ramp(j+wh)
            vert=ramp(j-wh):ramp(j+wh)
            p[i,j] = averaging!(view(v,vert,horz), view(ind,vert,horz))
        end
    end
end




# 時間平均するチャネル周波数応答(CFR)の推定
function get_ave_cfr(params, rx_frame, tx_frame; w_size=size(rx_frame,2))
    nblock = size(rx_frame,2) # OFDMシンボル数
    pind = params.pilot_indices # pilot indices
    div = Int(ceil(nblock/w_size)) # ブロックを分割
    cfr_temp = zeros(eltype(ComplexF64), params.fft_size, div, params.n_dim[2], params.n_dim[1]) # 平均チャネル応答配列
    period = size(params.pilot_order,1)*params.pilot_interval # パイロット送信周期(時間)
    timeinds = [Vector(1+params.pilot_interval*(i-1):period:nblock) for i in axes(params.pilot_order,1)]
    for m in axes(cfr_temp, 4) # Nr
        for n in axes(cfr_temp,3) # Nt
            for j in 1:div # フレーム分割
                for p in axes(params.pilot_order,1)# パイロット周期
                    l = params.pilot_order[p,n] # column of pilot indices!
                    blen = j*w_size < nblock ? j*w_size : nblock # 最大値の制限
                    tinds = filter(x-> (j-1)*w_size<x<=blen, timeinds[p])
                    if l > 0 # パイロットが送信されている場合
                        for i in view(pind,:,l)
                            cfr_temp[i,j,m,n] = sum(view(rx_frame,i,tinds,m) ./ view(tx_frame,i,tinds,n))/length(tinds) # average
                        end
                    end
                end
            end
        end
    end
    return cfr_temp
end






"""
推定チャネルベクトルと観測パイロット信号との誤差から雑音分散を推定する
"""
function noise_var_estimate(pilot::Pilot, Yp::AbstractArray{T,4}, Ypa::AbstractArray{T,3}) where T
    nvar = zeros(size(Yp,3))
    time_table = OFDM.make_time_table(pilot, size(Yp,2))
    Ntx = size(Yp,4)
    for i in axes(Yp,3) # Nrx
        total_sum = 0.0  # 二乗誤差
        M = 0   # 時間方向平均したサンプル数
        for j in axes(Yp,4) # Ntx
            for id in 1:pilot.n_group
                fidxs = pilot.gidxs[id] # idに対応する周波数index
                for k in fidxs # freq
                    N = 0
                    temp = 0.0
                    for l in axes(Yp,2) # time index
                        id != time_table[j,l] && continue
                        temp += abs2(Yp[k,l,i,j] - Ypa[k,i,j]) # 2乗誤差
                        N += 1
                    end
                    total_sum += temp / (N-1)
                    M += 1 # シンボル数
                end
            end
        end
        nvar[i] = total_sum / M # 不偏分散
    end
    nstats = length(get_time_idxs(pilot, size(Yp,2))) / pilot.n_group # 時間方向の統計量
    nvar ./= abs2(pilot.symbols[1]) * nstats # 推定分散
    return nvar
end
#=
function noise_var_estimate(pilot::Pilot, Yp::AbstractArray{T,4}, Ypa::AbstractArray{T,3}) where T
    nvar = zeros(size(Yp,3))
    time_table = OFDM.make_time_table(pilot, size(Yp,2))
    Ntx = size(Yp,4)
    for i in axes(Yp,3) # Nrx
        total_sum = 0.0  # 二乗誤差
        M = 0   # 時間方向平均したサンプル数
        for j in axes(Yp,4) # Ntx
            for k in axes(Yp,2) # time
                id = time_table[j,k]
                id == 0 && continue
                fidxs = pilot.gidxs[id] # idに対応する周波数index
                for l in fidxs  # freq
                    total_sum += abs2(Yp[l,k,i,j] - Ypa[l,i,j]) # 2乗誤差
                    M += 1
                end
            end
        end
        nvar[i] = total_sum / (M-1) # 不偏分散
    end
    nstats = length(get_time_idxs(pilot, size(Yp,2))) / pilot.n_group # 時間方向の統計量
    nvar ./= abs2(pilot.symbols[1]) * nstats # 推定分散
    return nvar
end
=#
