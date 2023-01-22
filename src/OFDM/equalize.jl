# equalizer.jl

mutable struct Time end
mutable struct Freq end

struct Equalizer{domain}
end

function Equalizer(domain=:Time)
    if domain == :Freq
        Equalizer{Freq}()
    elseif domain == :Time
        Equalizer{Time}()
    else
        error()
    end
end







"""
    (::Equalizer{Time})(ofdm, y, h, method)

Time-domain Equalization

====
    Arguments:
        ofdm   : OfdmType
        rxsig  : received signal in time domain
        h      : channel response vector
        method : detection scheme
====
    Returns:
        rxdata : equalized data

"""
function (::Equalizer{Time})(ofdm::Ofdm, rx_sig, h, nvar; estimator=:ZF)
    @unpack n_fft, n_gi, nTs, ndim, gitype, frame = ofdm
    ntx, nrx = ndim

    # Stack the received vector and Make Channel Toeplitz Matrix
    if ofdm.gitype isa CP
        y = stack(rx_sig, n_fft, n_gi, nrx)
        y = @view y[n_gi*nrx+1:end, :]
        H = make_channel_matrix(h, n_fft, n_fft, n_gi, circular=(gitype isa CP))
    else
        y = stack(rx_sig, n_fft, n_gi, nrx)
        H = make_channel_matrix(h, n_fft+n_gi, n_fft, n_gi, circular=(gitype isa CP))
    end

    # H * F
    #=
    for i in 1:nrx
        @views ifft!(H[i:nrx:size(H,1),:], 2)
    end
    H *= sqrt(n_fft*ntx) # scaling
    =#

    # Calculate a Weight Matrix
    if estimator==:MF
        W = H
    elseif estimator==:ZF
        W = H*inv(H'H)
    elseif estimator==:MMSE
        W = H*inv(H'H + nvar*I)
    end

    # Time Domain Equalization
    xhat = W'y

    # Convert into Frequency domain
    for i = 1:ntx
        @views fft!(xhat[i:ntx:end,:], 1) ./ sqrt(n_fft)
    end

    # Output
    xhat = reshape(xhat, ntx, n_fft, :)
    didxs = get_idxs(frame, :data); ndata = frame.on[:data]
    rx_data = @views xhat[:, didxs] |> vec
    return @views rx_data[1:ndata]
end



"""
    equalize(Y, H, N0=0.0; method=:ZF, output=:hard, option=nothing)
    equalize(ofdm::Ofdm, rx_frame::AbstractArray{T,3}, CFR; modt, N0=nothing, method=:ZF, required=:hard) where T

OFDMシンボル等化
"""
function equalize(ofdm::Ofdm, rx_frame::AbstractArray{T,3}, CFR::AbstractArray{T,3}, N0=nothing; at=axes(rx_frame,2), output_type=:hard, kwargs...) where T
    # フレームデータ
    ntx, nrx, nsym = ofdm.ndim..., size(rx_frame, 2)

    # 受信データ列
    if length(at) <= 1
        data_idxs = get_idxs_list(ofdm.frame, :data)[at] # データシンボルindex
        Y = @view rx_frame[data_idxs, at, :]
        H = @view CFR[data_idxs,:,:]
    else
        data_idxs = get_idxs(ofdm.frame, :data) # データシンボルindex
        Y = @view rx_frame[data_idxs, :]
        CFR = reshape(CFR, size(CFR,1), 1, size(CFR,2), size(CFR,3))
        CFR = repeat(CFR, 1, nsym, 1, 1)
        H = @view CFR[data_idxs,:,:]
    end

    # 雑音電力密度N0列(スカラorベクトル)
    if !isnothing(N0)
        ndims(N0)==2 && (N0 = @view N0[data_idxs])
        ndims(N0)==3 && (N0 = @view N0[:,data_idxs])
    end

    # signal_detectionを呼び出す
    outputs = signal_detection(Y, H, N0; dims=1, output_type=output_type, kwargs...)

    # データ数を検算して有効データのみ返す
    # numofdata = ofdm.frame.on[:data]
    if output_type == :hard
        rx_data = outputs |> vec
        return rx_data
    else
        rx_data = outputs[1] |> vec
        nvars   = outputs[2] |> vec
        return @views rx_data, nvars
    end
end


function equalize(Y, H, N0=nothing; idxs=Int[], mod_type, detector=ZF(), xhat=nothing, output_type=:hard)
    if !isempty(idxs)
        Y = @view Y[idxs,:]
        H = @view H[idxs, :, :]
    end
    if ndims(N0)==2
        N0 = @view N0[idxs,:]
    elseif ndims(N0)==3
        N0 = @view N0[idxs,:,:]
    end
    if !isnothing(xhat)
        xhat = @view xhat[idxs,:]
    end
    return signal_detection(Y, H, N0; dims=1, mod_type=mod_type, detector=detector, xhat=xhat, output_type=output_type)
end
