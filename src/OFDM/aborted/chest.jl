
"""
パイロットチャネル推定
"""
function chest(ofdm::Ofdm, rx_frame::AbstractArray{T}, tx_frame::AbstractArray{T}
                ;average=false, interpolation=true, method=:linear) where T

    @assert size(rx_frame,2) == size(tx_frame,2)
    isnothing(ofdm.pilot) && error("pilotは割り当てられていないためチャネル推定できません。")
    ntx, nrx = ofdm.ndim
    pilot = ofdm.pilot
    alist = ofdm.frameinfo.alist
    nfft, nsym = size(rx_frame)[1:2]

    # CFR配列確保
    CFRs = zeros(T, nfft, nsym, nrx, ntx)
    _pilot_divide!(CFRs, rx_frame, tx_frame, alist)
    interpolation && interpolate!(CFRs, method)

    if average
        return _average!(CFRs, pilot)
    else
        return CFRs
    end
end


# パイロット除算
function _pilot_divide!(CFRs, rx_frame, tx_frame, alist)
    for j in axes(CFRs,2)
        for i in axes(CFRs,1)
            alist[i,j] <= 0 && continue
            p = alist[i,j]
            for q in axes(CFRs,3)
                CFRs[i,j,q,p] = rx_frame[i,j,q] / tx_frame[i,j,p]
            end
        end
    end
end

struct Time end
struct Freq end
struct TimeFreq end

# パイロット補間
function interpolate!(pilot::Pilot, CFRs, dir=:freq)
    ord = pilot.order
    dinds = ofdm.dinds

    if dir==:freq
        _interpolate(Freq(), )
    elseif dir == :time
        _interpolate(Time(), )
    elseif dir == :timefreq
        _interpolate(Freq(), method=method)
        _interpolate(Time(), method=:identity)
    end
end

# 周波数方向へ補間
function _interpolate!(::Freq, pilot::Pilot, CFRs)
    ord = pilot.order
    stop = size(CFRs, 2)
    for q in axes(CFRs,3)
        for p in axes(CFRs,4)
            for id in 1:n_group
                start = @views t0 + findfirst(x->x==id, ord[p,:])
                knots = start:step:stop
                for i in pilot.idxs[id]
                    @views lininterp!(CFRs[i,:,p,q], knots)
                end
            end
        end
    end
end






"""
パイロットチャネル推定結果のブロック平均をとる
"""
function average(CFRs::AbstractArray, pilot::Pilot, n_blocks)
    n_block < 2 && @error "filter_sizeは2以上である必要があります"
    n_frame = size(CFRs,2) # フレームサイズ
    divs = div(n_frames, n_blocks) # フレーム分割数
    nfft, n_frames, n_rx, n_tx = size(CFRs)
    ave_cfrs = zeros(eltype(CFRs), nfft, divs, n_rx, n_tx) # 出力配列
    time_inds = gen_pilot_time_ids(pilot, n_frame)
    pilot.n_cluster==0 && @error "パイロットが挿入されていません。"
    for q in axes(CFRs,3)
        for p in axes(CFRs,4)
            for id in 1:pilot.n_cluster
                m = 1; n = 1; count = 0
                pinds = pilot.indices[id]
                while n<=n_frames
                    if id == time_inds[n,p]
                        ave_cfrs[pinds,m,q,p] .+= CFRs[pinds,n,q,p] # add
                        count += 1
                    end
                    if n%n_block==0 || n==n_frames # 平均を取って次ブロックに移動
                        count==0 && @warn "n_blocksが不適切です。"
                        ave_cfrs[pinds,m,q,p] ./= count # average
                        count = 0
                        m += 1
                    end
                    n += 1
                end
            end
        end
    end
    ave_cfrs
end


"""
チャネル推定
"""
function channel_estimate!(rxframe::AbstractMatrix{T}, txframe::AbstractMatrix{T}, weight) where T
    for q in axes(rxframe,3)
        for p in axes(txframe,3)
            _channel_update!(Ĥ, Y, X, pinds, weight)
        end
    end
end

# チャネル係数更新
function _channel_update!(Ĥ, Y, X, pinds, weight)
    for k in pilot_inds
        Ĥ_new = Y[k] ./ X[k]
        Ĥ[k] = (1-weight)*Ĥ_new + weight*Ĥ[k]
    end
end
