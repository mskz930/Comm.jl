

struct Slicer
  modulator
end

#=
(s::Slicer)(inp::AbstractVector{T}) where {T}
  out = zeros(T, length(inp))
  η = normfactor(s.modulator)
  for n = 1:length(inp)
    out[n] = slice(modulator, inp[n]) / η
  end
end
=#

# Slicer(qam::Qam) = Slicer{typeof(qam)}(normfactor(qam))

# (s::Slicer{Qam{2}})(r) = slice_2qam.(r)
# (s::Slicer{Qam{4}})(r) = (r *= s.η; slice_4qam(r)/s.η)
# (s::Slicer{Qam{16}})(r) = (r *=s.η; slice_16qam(r)/s.η)




#=
slice(pam::Pam, rx_data) = slice(PAM, rx_data, M=pam.M, isnorm=pam.isnorm)


# QAM slicer
function slice(::Type{QAM}, rx_data::T, M, isnorm) where T<:Number
    if M == 2
        sliced_data = slice_qam2(rx_data)
    elseif M == 4
        isnorm && (rx_data *= sqrt(2))
        sliced_data = slice_qam4(rx_data)
        isnorm && (sliced_data /= sqrt(2))
    elseif M == 16
        isnorm && (rx_data *= sqrt(10))
        sliced_data = slice_qam16(rx_data)
        isnorm && (sliced_data /= sqrt(10))
    elseif M == 64
        isnorm && (rx_data *= sqrt(41))
        sclies_data = slice_qam64(rx_data)
        isnorm && (sliced_data /= sqrt(41))
    else
        error("M=", M, "の値が適切でありません.")
    end
    sliced_data
end

function slice(::Type{QAM}, rx_data::T, M, isnorm) where T<:AbstractArray
    if M == 2
        sliced_data = slice_qam2.(rx_data)
    elseif M == 4
        isnorm && (rx_data .*= sqrt(2))
        sliced_data = slice_qam4.(rx_data)
        isnorm && (sliced_data ./= sqrt(2))
    elseif M == 16
        isnorm && (rx_data .*= sqrt(10))
        slices_data = slice_qam16.(rx_data)
        isnorm && (sliced_data ./= sqrt(10))
    elseif M == 64
        isnorm && (rx_data .*= sqrt(41))
        sclies_data = slice_qam64.(rx_data)
        isnorm && (sliced_data ./= sqrt(41))
    else
        error("M=", M, "の値が適切でありません.")
    end
    sliced_data
end


# Qamオブイェクトからslice(QAM)を呼び出すラッパー
slice(qam::Qam, rx_data) = slice(QAM, rx_data, M, qam.isnorm)
=#
