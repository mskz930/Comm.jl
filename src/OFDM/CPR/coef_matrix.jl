
make_coef_mat(ofdm::Ofdm, h; domain=:freq) = make_coef_mat(h, ofdm.n_fft, ofdm.n_gi; domain)

# 干渉行列の生成
function make_coef_mat(h, n_fft, n_gi; domain=:time)
    Hici = zeros(eltype(h), n_fft, n_fft, size(h,2), size(h,3))
    Hisi = zeros(eltype(h), n_fft, n_fft, size(h,2), size(h,3))
    if (size(h,1)-1)>n_gi
        make_coef_mat!(Hisi, Hici, h, n_fft, n_gi, domain)
    end
    Hisi, Hici
end

# 干渉行列の生成
function make_coef_mat!(Hisi, Hici, h, n_fft, n_gi, domain=:time)
    for n in axes(Hisi,4)
        for m in axes(Hisi,3)
            h_, Hisi_, Hici_ = @views h[:,m,n], Hisi[:,:,m,n], Hici[:,:,m,n]
            _make_coef_mat_in_time!(Hisi_, Hici_, h_, n_gi)
            if domain==:freq
                fft!(Hisi_, 1); Hisi ./= sqrt(n_fft)
                ifft!(Hisi_, 2); Hisi .*= sqrt(n_fft)
                fft!(Hici_, 1); Hici ./= sqrt(n_fft)
                ifft!(Hici_, 2); Hici .*= sqrt(n_fft)
            end
        end
    end
end

# テプリッツ行列(干渉成分)の生成
function _make_coef_mat_in_time!(Hisi, Hici, h, G)
    N, L = size(Hisi,1), (size(h,1)-1) # FFT数, 遅延タップ数

    # 時間領域Hisiの作成
    for i in 1:L-G
        for j in i:L-G
            Hisi[i, end-j+i] = h[1+G+j]
        end
    end
    # 時間領域Hiciの作成
    for i in 1:L-G
        for j in i+G:L
            Hici[i, end-j+i] = h[1+j]
        end
    end
end



# ISI係数行列(周波数領域)の計算:
function _make_coef_mat_in_freq!(Hi, h, W, G)
    N, L = size(Hi,1), size(h,1) - 1
    temp = complex(0.0)
    d = 0
    for k in 0:N-1
        for m in 0:N-1
            temp = complex(0.0)
            d = (k-m)
            if d == 0
                for l in G+1:L
                    Hi[k+1,m+1] += (l-G) * h[l+1] * W[m+1, l+1] / N
                end
            else
                for l in G+1:L
                    if d > 0
                        for n in 0:l-G-1
                            temp += W[n+1,d+1]
                        end
                    else
                        for n in 0:l-G-1
                            temp += conj(W[n+1,-d+1])
                        end
                    end
                    Hi[k+1,m+1] += temp * h[l+1] * W[m+1, l+1] / N # exp(-2π*m*l/N)
                end
            end
        end
    end
end

# ISI係数(周波数領域)計算
function isi_coef_calc(h::AbstractVector, F::AbstractMatrix, L, G, N, k, m)
    d = k-m
    coef = zero(eltype(h))
    if d == 0
        for l in G+1:L
            iszero(h[l+1]) && continue
            coef += (l-G) * h[l+1] * F[m+1,l+1]
        end
        return coef * conj(F[m+1,G+1]) / N
    else
        for l = G+1:L
            iszero(h[l+1]) && continue
            temp = zero(eltype(h))
            if d>0
                for n = 0:l-(G+1)
                    temp += F[n+1,d+1] # sum_{n=0}^{l-G-1} exp(-2πn(k-m)/N)
                end
            else
                for n = 0:l-(G+1)
                    temp += conj(F[n+1,-d+1]) # sum_{n=0}^{l-G-1} exp(-2πn(k-m)/N)
                end
            end
            coef += temp * h[l+1] * F[m+1,l+1] # 1/N * h(l) * exp(-2πml/N)
        end
        return coef * conj(F[m+1,G+1]) / N
    end
end

#ICI係数(周波数)の計算
function ici_coef_calc()
end


function ici_channel_matrix(H::AbstractMatrix{T}, G, L) where T
    H0 = copy(H)
    H1 = zero(H)
    R = (L-1)-G
    for i in 1:R
        for j in 0:R-i
            H1[i, end-(L-1)+i+j] = H[i,end-(L-1)+i+j]
            H0[i, end-(L-1)+i+j] = zero(T)
        end
    end
    H0, H1
end

#=
"""
ICI, ISI行列の生成(非ゼロな要素に対してのみ素直に計算)
"""
function make_coef_matrix_in_time(h::AbstractArray, N::Integer, G::Integer, W::AbstractMatrix)
    # h: channel impluse response vector(possible multidimension)
    # N, G, W: FFT点数, GI長, DFT行列

    # 事前計算
    n_tx = size(h,3); n_rx = size(h,2);
    L = size(h,1)
    H_ici = zeros(Complex, N, N, size(h,2), size(h,3)) # ICI係数行列
    H_isi = zeros(Complex, N, N, size(h,2), size(h,3)) # ISI係数行列
    ls = findall(x -> x > 0, abs.(view(h,:,1,1))) # channel coef indices
    ls = filter(x -> x > G+1, pind) .- 1

    temp = complex(0.0)

    for q in 1:n_rx
        for k in 0:N-1
            for m in 0:N-1
                if k != m # k!=m
                    for l in ls # nonzero index over G
                        temp = zero(Complex)
                        for n in 0:(l-G)-1
                            temp += (k-m) > 0 ? W[n+1,(k-m)+1] : conj(W[n+1,(m-k)+1])
                        end
                        for p in 1:n_tx
                            H_ici[k+1,m+1,q,p] += temp * h[l+1,q,p] * W[m+1,l+1]/N
                            H_isi[k+1,m+1,q,p] += temp * h[l+1,q,p] * W[m+1,l-G+1]/N
                        end
                    end
                else # k==m
                    for p in 1:n_tx
                        temp = zero(Complex)
                        for l in ls
                            temp += (l-G) * h[l+1,q,p] * W[k+1,l+1] # (l-G) * h[l,q,p] * exp(-j2πkl/N)
                        end
                        H_ici[k+1,m+1,q,p] = temp/N
                        H_isi[k+1,m+1,q,p] = temp * conj(W[k+1,G+1])/N
                    end
                end
            end
        end
    end
    H_ici, H_isi
end
=#
