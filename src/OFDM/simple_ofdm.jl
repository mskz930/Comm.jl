# 
using Random
using FFTW

using PyPlot; pygui(false)
N_FFT = 64
N_SYM = 10
N_GI = 16
N_TX = 2

bps = 2
n_bits = N_FFT * N_SYM * N_TX * bps
bits = bitrand(n_bits) |> Vector
symbols = map(x -> x > 0 ? 1 : -1, bits)
symbols = reshape(symbols, bps, :)
symbols = symbols[1,:] .+ im * symbols[2,:]

struct OFDM
    
end

function ofdmmod(symbols, n_fft, n_gi, n_tx=1)
    # to tensor
    ofdm_symbols = reshape(symbols, n_tx, n_fft, :)
    # IFFT
    ofdm_signal = ifft(ofdm_symbols, 2) / sqrt(n_fft)
    # add GI
    ofdm_signal = [ofdm_signal[:, end-n_gi+1:end, :]; ofdm_signal]
    return vec(ofdm_signal)
end


function odfmdemod()
end



plt.figure()
plt.plot(real(signal))
plt.plot(imag(signal))
display(gcf())