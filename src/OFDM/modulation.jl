# ofdmmod.jl

export ofdm_mod, ofdm_demod


# OFDM変調
function mod(ofdm::OfdmMod, tx_frame)
    @unpack n_fft, n_gi, n_dim, gitype, giloc = ofdm
    ifft_out = sqrt(ofdm.n_fft) .* ifft(tx_frame, 1) # N-point IFFT
    ofdm_sig = _giadd(ifft_out, n_gi, gitype, giloc) # ガードインターバル付与
    return reshape(ofdm_sig, :, n_dim[1])
end

ofdm_mod = mod # alias


#=
GI(Prefix/Postfix)を付加する
    CP => Cyclic Prefix
    ZP => Zero Postfix
=#
_giadd(x, G, ::CP, ::Prefix) = @views vcat(x[end-G+1:end,:,:], x)
_giadd(x, G, ::ZP, ::Prefix) = vcat(zeros(eltype(x), G, size(x,2), size(x,3)), x)

_giadd(x, G, ::CP, ::Postfix) = @views vcat(x, x[end-G+1:end, :, :])
_giadd(x, G, ::ZP, ::Postfix) = vcat(x, zeros(eltype(x), G, size(x,2), size(x,3)))


# OFDM復調
function demod(ofdm::OfdmMod, rx_sig::AbstractArray{T}) where T
    nrx = ofdm.n_dim[2]
    n_Ts = symlen(ofdm)
    rx_sig = reshape(rx_sig, n_Ts, :, nrx)
    girem = rx_sig[ofdm.n_gi+1:ofdm.n_fft+ofdm.n_gi, :, :]
    rx_frame  = 1/sqrt(ofdm.n_fft) .* fft!(girem, 1)
    return rx_frame
end

ofdm_demod = demod

chest(ofdm::OfdmMod, rx_frame, tx_frame) = OFDM.Chest.channel_estimate(rx_frame, tx_frame; ofdm.chest...)

