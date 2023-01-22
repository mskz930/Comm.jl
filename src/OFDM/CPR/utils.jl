
"""
    channel_compensation(cfr, h, n_fft, n_gi)
周波数応答を補正する
"""
function channel_compensation(cfr, h, N, G)
    H = copy(cfr)
    H = channel_compensation!(H, h, N, G)
    H
end

function channel_compensation!(H, h, N, G)
    L = size(h,1)-1
    T = ComplexF64
    tmp = zero(T)
    for q in axes(H,2)
        for p in axes(H,3)
            for k in 0:N-1 # freq
                tmp = zero(T)
                # チャネル係数の補正項
                for l in G+1:L
                    tmp += (l-G)/N * h[l+1,q,p] * exp(-im*(2π*k*l/N))
                end
                H[k+1,q,p] -= tmp
            end
        end
    end
    H
end


"""
    calc_sinr(ofdm, h, σn2=0.01; scale)
SINRとそのときの等価雑音分散を計算する

arguments:
    N  : FFT
    G  : Guard Interval
    Nd : Data Subcarriers
    N0 : Noise Power Density
    scale : :linear, :log10

returns:
    SINRs: SINR[dB]
    nvars: 等価雑音電力

"""
function calc_SINR(h::AbstractArray, N::Int, G::Int, Nd::Int, nvar::Real; scale=:linear)
    svar = Nd/N
    Nr, Nt = size(h,2), size(h,3)
    L = size(h,1)

    SINR = zeros(Float64, Nr); invar = zeros(Float64, Nr)
    for i in 1:Nr
        Ps = 0.0 # 平均信号電力
        Pin = N*nvar # 平均干渉電力
        for j in 1:Nt
            for l in 1:L
                if l < G+2
                    Ps += svar*abs2(h[l,i,j]) * N
                else
                    Ps += svar*abs2(h[l,i,j]) * (N-l+G+1)
                    Pin += svar*abs2(h[l,i,j]) * (l-G-1) # 平均干渉電力=平均受信電力*サンプル数
                end
            end
        end
        SINR[i] = scale==:linear ? Ps/Pin : 10log10(Ps/Pin)
        invar[i] = Pin/Ps
    end
    SINR, invar
end

calc_SINR(ofdm::Ofdm, h::AbstractArray, N0::Real=0.0; scale=:linear) = calc_SINR(h, ofdm.n_fft, ofdm.n_gi, length(ofdm.idxs[:data]), N0; scale)
