
const CF64 = ComplexF64

# get_ici_matrix: ICI行列取得
function get_coef_matrix(h::AbstractArray, N::Integer, G::Integer, W::AbstractMatrix)
    # h: channel impluse response vector(possible multidimension)
    # N, G, W: FFT点数, GI長, DFT行列

    # 事前計算
    n_tx = size(h,3); n_rx = size(h,2);
    L = size(h,1)
    H_ici = zeros(CF64, N, N, size(h,2), size(h,3)) # ICI係数行列
    H_isi = zeros(CF64, N, N, size(h,2), size(h,3)) # ISI係数行列
    ls = findall(x -> x > 0, abs.(view(h,:,1,1))) # channel coef indices
    ls = filter(x -> x > G+1, pind) .- 1

    # 係数行列を求める計算
    for q in 1:n_rx
        for k in 0:N-1
            for m in 0:N-1
                if k != m # k!=m
                    for l in ls
                        temp = zero(CF64)
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
                        temp = zero(CF64)
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

function get_coef_matrix(cir::AbstractArray, n_fft::Integer, n_gi::Integer, domain::Symbol=:time)
    # 損失チャネル行列の生成
    isi_mtx = zeros(eltype(cir),n_fft,n_fft,size(cir,2), size(cir,3)) # isi行列
    L = size(cir,1) # チャネル長(パス数)
    Q = L-n_gi-1
    T = view(isi_mtx, 1:Q, n_fft-Q+1:n_fft, :, :)
    for q in axes(cir,2)
        for p in axes(cir,3)
            for m in 1:Q
                for n in m:Q
                    l = n-m
                    T[m,n,q,p] = cir[end-l,q,p] #
                end
            end
        end
    end
    ici_mtx = circshift(isi_mtx, (0,-n_gi))
    if domain==:freq
        fft!(isi_mtx, 1); fft!(ici_mtx, 1); # 列方向でFFT
        ifft!(isi_mtx, 2); ifft!(ici_mtx, 2); # 行方向でIFFT
    end
    ici_mtx./n_fft, isi_mtx./n_fft
end


function coef_matrix_calc(Hici, Hisi, h::AbstractVector, n_fft, n_gi, F)
    for k in 0:n_fft-1
        for m in 0:n_fft-1
            d = k-m
            if k !== m
                for l in pind
                    b = zero(CF64)
                    for n in 0:(l-n_gi)-1
                        b += d > 0 ? F[n+1,d+1] : conj(F[n+1,-d+1])
                        # b += exp(-im*(2*pi*d*n)/n_fft)
                    end
                    for p in 1:n_tx
                        Hici[k+1,m+1,q,p] += b * h[l+1,q,p] * F[m+1,l+1] / n_fft
                        Hisi[k+1,m+1,q,p] += b * h[l+1,q,p] * F[m+1,l-n_gi+1] / n_fft
                    end
                end
            else

                    b = zero(CF64)
                    for l in pind # 遅延軸
                        b += (l-n_gi)/n_fft* h[l+1,q,p] * F[k+1,l+1]
                    end
                    ici_mtx[k+1,m+1,q,p] = b
                    isi_mtx[k+1,m+1,q,p] = b*conj(F[k+1,n_gi+1])
                end
            end
        end
    end
end
