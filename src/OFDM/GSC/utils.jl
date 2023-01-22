# GSC/utils.jl

# 受信信号ベクトルrのスタック
function _stackvec!(v1, v2)
    n = size(v2,1)
    for k in axes(v2, 3)
        for j in axes(v2, 2)
            for i in axes(v2, 1)
                v1[i+n*(k-1),j] = v2[i,j,k]
            end
        end
    end
end

# stackしたvector信号を多次元配列に戻す
function _destackvec!(v1, v2)
    ntx, nfft, _ = size(v1)
    for k in axes(v1,3) # time
        for j in axes(v1,2) # freq
            for i in axes(v1,1) # space
                v1[i,j,k] = v2[j+(i-1)*nfft,k]
            end
        end
    end
end


# 行列のstack
function _stackmat(H)
    N, Nr, Nt = size(H,1), size(H,3), size(H,4)
    if Nt == 1
        Hstack = reshape(H, :, Nr*Nt)
        Hstack = transpose(Hstack)
        Hstack = reshape(Hstack, Nr*N, Nt*N)
    else
        Hstack = permutedims(H, [3,1,4,2])
        Hstack = reshape(Hstack, Nr*N, Nt*N)
    end
    return Hstack
end
function _stackmat!(H1, H2)
    nfft, nrx, ntx = size(H2)[2:4]
    for l in axes(H2,4)
        for k in axes(H2,3)
            for j in axes(H2,2)
                for i in axes(H2,1)
                    H1[i+nfft*(k-1),j+nfft*(l-1)] = H2[i,j,k,l]
                end
            end
        end
    end
end
