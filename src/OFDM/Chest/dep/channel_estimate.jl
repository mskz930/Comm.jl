# channel_est.jl


"""
パイロットシンボルによるOFDMチャネル推定
"""
function channel_estimate(pars::OfdmMod, rx_frame::AbstractArray{T}, tx_frame::AbstractArray{T};
                     interp=false, method=:linear, average=true, domain=:time) where T <: Number
    # FFT数, ブロック数, 受信アンテナ, 送信アンテナ
    n_fft, n_blocks, n_rx, n_tx = (size(rx_frame), size(tx_frame,3))
    # CFR推定結果配列確保
    cfrs = zeros(T, N, n_blocks, n_rx, n_tx)
    #
    @views divide!(cfrs, rx_frame, tx_frame, inds) # パイロット除算
    interp && interpolate!(cfrs, method) # 補間
    average && average!(cfrs, )
    return cfr_est
end



function channel_estimation(params, rx_frame, tx_frame; interp=true, average=false, w_size=size(rx_frame,2))
    # FFT数, ブロック数, 受信アンテナ, 送信アンテナ
    nfft, nblocks, Nr, Nt = (size(rx_frame), size(tx_frame,3))
    # CFR推定用配列
    cfr = zeros(eltype(rx_frame), N, nblocks, Nr, Nt)
    interval = params.pilot.interval # 挿入間隔(時間)
    pattern = size(params.pilot.order,1) # パイロット送信パターン数
    period = pattern*interval # パイロット挿入周期
    pilot_indices = params.pilot.indices # パイロットindex

    for n in axes(cfr,4) # Nt
        for m in axes(cfr,3) # Nr
            for p in 1:pattern # パイロット周期
                l = params.pilot.order[p,n] # column of pilot indices!
                if l > 0 # パイロットが送信されている場合
                    knots = Vector(1+(p-1)*interval:period:nblocks) # 時間index
                    for i in pilot.indices[l]
                        for j in knots
                            cfr[i,j,m,n] = rx_frame[i,j,m] / tx_frame[i,j,n] # パイロット除算
                        end
                        if average # 時間平均化
                            @views cfr[i,:,m,n] = pilot_averaging(cfr[i,:,m,n], knots, w_size)
                        elseif interp # 補間
                            lininterp!(view(cfr,i,:,m,n), knots) # 線形補間(時間方向)
                        end
                    end
                end
            end
            if interp
                pinds = sort!(union(pilot.indices...)) # パイロットindex集合
                for j in axes(cfr,2)
                    lininterp!(view(cfr,:,j,m,n), pinds) # 線形補間(周波数補間)
                end
            end
        end
    end
    return cfr
end

# パイロット平均化
function average(input::AbstractVector, inds, w_size=1)
    @assert w_size%2!==0
    if length(inds)>1
        output = Array{eltype(input)}(undef,length(input))
        w_size_half = Int((w_size-1)/2)
        temp_inds = filter(x -> 1-w_size_half<=x<=1+w_size_half, inds)
        mean_val = mean(view(input,temp_inds)) # 初期平均値
        ninds = temp_inds[end] < inds[end] ? inds[length(temp_inds)+1] : length(input)+w_size_half
        for i in 1:length(input)
            if i <= inds[end]+w_size_half
                if temp_inds[1]<=i-w_size_half || ninds <= i+w_size_half
                    temp_inds = filter(x -> i-w_size_half<=x<=i+w_size_half, inds)
                    ninds = temp_inds[end] < inds[end] ? inds[length(temp_inds)+1] : length(input)+w_size_half
                    mean_val = mean(view(input,temp_inds))
                end
                output[i] = mean_val
            else
                output[i] = output[i-1]
            end
        end
        return output
    else
        return input
    end
end




# 逐次チャネル推定
function channel_estimation(params::Ofdm, rx_frame, tx_frame, current_cfr, j, frameinfo)
    if j%params.pilot.interval == 0
        cfr = zeros(eltype(rx_frame), size(rx_frame,1), size(rx_frame,2), params.Nr, params.Nt) # cfr推定
        interval = params.pilot.interval # 挿入間隔(時間)
        pattern = size(params.pilot.order,1)
        period = pattern*interval # パイロット挿入周期
        nblock = size(rx_frame,2) # ブロック長
        pilot_indices = params.pilot.indices # パイロットindex
        for n in axes(cfr,4) # Nt
            for m in axes(cfr,3) # Nr
                for p in 1:pattern # パイロット周期
                    l = params.pilot.order[p,n] # column of pilot indices!
                    if l > 0 # パイロットが送信されている場合
                        for i in view(pilot.indices,:,l)
                            for j in 1+interval*(p-1):period:nblock
                                cfr[i,j,m,n] = rx_frame[i,j,m] / tx_frame[i,j,n] # パイロット除算
                            end
                            knots = collect(p:period*params.pilot.interval:nblock) # パイロットが含まれる時間インデックス
                            lininterp!(view(cfr,i,:,m,n), collect(knots)) # 線形補間(時間方向)
                        end
                        knots = vec(pilot.indices')
                        for j in axes(cfr,2)
                            lininterp!(view(cfr,:,j,m,n), knots) # 線形補間(周波数補間)
                        end
                    end
                end
            end
        end
    else
    end
    return cfr
end


# パイロット除算
function divide!(Ĥ, Y, X, frameinfo)
    for n in axes(Ĥ,4) # Nt
        for m in axes(Ĥ,3) # Nr
            for j in axes(Y,2)
                for i in axes(Y,1)
                    if frameinfo[i,j,n] > 0
                        Ĥ[i,j,m,n] = Y[i,j,m] / X[i,j,n] # パイロット除算
                    end
                end
            end
        end
    end
end

# パイロット除算
# inds : (freq,time) indices
function divide!(Ĥ, Y, X, inds)
    n_rx, n_tx = size(Ĥ)[3:4]
    for (k,l) in Tuple.(inds)
        for n in 1:n_tx
            for m in 1:n_rx
                Ĥ[k,j,m,n] = Y[k,l,m] / X[k,l,n]
            end
        end
    end
end
