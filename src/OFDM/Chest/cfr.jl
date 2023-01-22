# cfr.jl

"""
    to_cfr(cir::AbstractArray{T}, fft_size::Integer; padding=true, n_repeats=1) where T

CIRベクトル列をCFRベクトル列に変換する

    Arguments:
        ofdm : Ofdm Type
        cir  : Channel Impulse Response
             時不変チャネルの場合: (L, Nrx, Ntx)
             時変動チャネルの場合: (L, Nsym, Nrx, Ntx)
    Returns:
        cfr  : Channel Frequency Response
"""
function to_cfr(ofdm::Ofdm, cir; padding=true)
    @unpack n_fft, n_gi = ofdm
    n_fft < size(cir,1) && @error "fft_sizeはチャネル長Lより大きい必要があります."
    if n_dims(cir) == 4
        # 時変動チャネルの場合時間方向平均をとる
        cir = _time_average(cir, n_fft+n_gi)
    end
    padding && (cir = zero_padding(cir, n_fft-size(cir, 1), dims=1)) # zero-pading
    if eltype(cir) <: Real
        cir = convert(Array{ComplexF64}, cir)
    end
    cfr = fft!(cir, 1) # fft
    cfr
end

function to_cfr!(cfr, cir)
    for p in axes(cir,2)
        for q in axes(cir,2)
            cfr[1:size(cir,1),q,p] .= cir[:,q,p]
            cfr[size(cir,1)+1:end,q,p] .= zero(eltype(cfr))
            @views fft!(cfr[:,q,p])
        end
    end
    cfr
end

# 時間平均
function _time_average(arr, per_sample)
    avg_size = cld(size(arr,2), per_sample)
    avg_out  = zeros(eltype(arr), size(arr,1), avg_size, size(arr)[3:end]...)
    for i in axes(avg_out,2)
        st = 1+(i-1)*per_sample
        en = i < size(avg_out,2) ? i*per_sample : size(arr,2)
        avg_out[:,i,:,:] = squeeze(mean(arr[:,st:en,:,:], dims=2))
    end
    avg_out
end


# CFRを時間平均する
function time_averaged_cfr(cfr::AbstractArray{T}, pind, frameinfo) where T
    ave_cfr = zeros(T, size(cfr,1), size(cfr,3), size(cfr,4)) # 平均化周波数応答
    for n in axes(cfr,4)
        for m in axes(cfr,3)
            for p in 1:size(pind,2)
                for i in view(pind,:,p)
                    ind = view(frameinfo,i,1:size(cfr,2),n)
                    knots = LinearIndices(ind)[ind .== p] # pilot index
                    if !isempty(knots)
                        ave_cfr[i,m,n] = _average(view(cfr,i,:,m,n), knots) # 線形補間
                    end
                end
            end
        end
    end
    ave_cfr
end
