# GSC.jl
module GSC

# packages
using LinearAlgebra
using FFTW
using SparseArrays

# modules
using ..OFDM
using ..Modulator
using ..OFDM.Chest: to_cfr
using ..CPR
using ...Comutils: dftmat
using ...Sigdet: columnsort
using ...Sigdet.MIMO

include("utils.jl")

"""
    equalize(ofdm::Ofdm, rx_frame::Frame, CIR, CFR, N0; method=:MMSE, approx=false)

GSC(Generalized Siderobe Canceller) 等化器:
    W_mmse = (I - B*(B'*Rin*B)^-1*B'*Rin)*D
"""
function equalize(ofdm::Ofdm, rx_frame, CIR, CFR,  N0; method=:MMSE, approx=false, option=())
    # パラメータ
    nfft, ngi, nTs, nsym = ofdm.nfft, ofdm.ngi, ofdm.nTs, size(rx_frame, 2)
    ntx, nrx = ofdm.ndim
    T = eltype(rx_frame)

    # 受信信号配列,　干渉行列の生成
    z = Array{T}(undef, nfft*nrx, nsym)
    _stackvec!(z, rx_frame)
    V, D, N0 = _make_mmse_filer(CIR, CFR, method, N0, approx, nfft, ngi, nrx, ntx)

    # 利得の計算
    if ntx == 1
        y = Array{T}(undef, ntx, nfft, nsym)
        D = _make_coherent_matrix(CFR)
        α = _gain(CFR)
        rxdata = (D'*V'*z) ./ α
        #D = _make_coherent_matrix(CFR)
        #_simo_detection(z, W, D, CFR)
        dinds  = get_all_idxs(ofdm, :data)
        return vec(rxdata[dinds])
    else
        y = Array{T}(undef, nrx, nfft, nsym)
        _destackvec!(y, V'z)
        dinds  = get_idxs_list(ofdm, :data)
        modt = get(option, :modt, nothing)
        rxdata = _mimo_detection(y, CFR, dinds, method, N0, modt)
        return vec(rxdata[:, get_all_idxs(ofdm, :data)])
    end
end

# GSCフィルタの生成
function _make_mmse_filer(CIR, CFR, method, N0, approx, nfft, ngi, nrx, ntx)
    Hisi  = _make_coef_matrix(CIR, nfft, ngi, nrx, ntx)
    D     = _make_coherent_matrix(CFR)
    B     = _make_block_matrix(CFR)
    Rin   = 2*Hisi*Hisi' + N0*I

    if approx
        F  = dftmtx(nfft, normalize=true)
        T  = _make_reduce_dim_mat(B, F, CIR, nfft, ntx, nrx)
        BT = B*T
        U  = inv(BT'*Rin*BT)*(BT'*Rin)
        V  = I - BT*U
    else
        U = inv(B'*Rin*B)*B'*Rin
        V  = I - B*U
    end

    if method==:MMSE || method==:SIC
        N0 = real(tr(V'*Rin)) / size(D,2)
        # N0 = real(tr(D'*V'*Rin*D)) / size(D,2) # MMSEの値
    end
    return V, D, N0
end

# 干渉行列(Hisi)の生成
function _make_coef_matrix(CIR, nfft, ngi, nrx, ntx)
    Hisi = Array{eltype(CIR)}(undef, nfft*nrx, nfft*ntx)
    Hisi_, _ = CPR.gen_coef_matrix(CIR, nfft, ngi, (ntx, nrx), domain=:freq)
    _stackmat!(Hisi, Hisi_)
    return Hisi
end

# 整合行列Dの生成
function _make_coherent_matrix(CFR)
    nfft, nrx, ntx  = size(CFR)
    D = spzeros(eltype(CFR), nfft*nrx, nfft*ntx)  # 整合フィルタ
    for j in axes(CFR,3)
        for i in axes(CFR,2)
            for k in axes(CFR,1)
                D[k+(i-1)*nfft,k+nfft*(j-1)] = CFR[k,i,j]
            end
        end
    end
    D
end

# ブロック行列Bの生成
function _make_block_matrix(CFR)
    nfft, nrx, ntx = size(CFR)
    nrem = nrx - ntx
    H = zeros(eltype(CFR), nrx, ntx)
    B = spzeros(eltype(CFR), nrx*nfft, nrem*nfft)
    bs = zeros(eltype(B), nrx, nrem)

    for k = 1:nfft
        # making vector b
        H .= @view CFR[k,:,:]
        # bs .= nullspace(H')
        if size(bs,2) == 1
            @views _make_b_vec!(bs, H')
        else
            @views _make_b_vecs!(bs, H')
        end

        # 代入
        for j in axes(bs,2)
            for i in axes(bs,1)
                B[k+(i-1)*nfft, k+(j-1)*nfft] = bs[i,j]
            end
        end
    end
    return B
end

# ブロックベクトルbの作成　
function _make_b_vec!(b, D)
    T = eltype(b)
    len = length(b)
    deg = size(D,2) - size(D,1) # 自由度

    # 自由度の分だけ係数を決定する
    _gauss_elimination!(D)
    b .= zero(T)
    b[end-deg+1:end] .= one(T)
    for i in length(b)-deg:-1:1
        for j in i+1:size(D,2)
            b[i] -= D[i,j] * b[j]
        end
    end
    b ./= norm(b) # normalize
end

# ブロックベクトルbの作成(複数ベクトルの場合　)
function _make_b_vecs!(bs, D)
    T = eltype(bs)
    len = size(bs,1)
    ndeg = size(D,2) - size(D,1)
    bs .= zero(T)

    _find_nullspace_basis!(D, bs)
    _orthonormalize_by_schmidt!(bs)
end

# 行列の右Null空間の基底を求める
function _find_nullspace_basis!(D, bs)
    m, n = size(D)
    T = eltype(D)

    _gauss_elimination!(D)

    for j in axes(bs,2)
        for i in 1:m
            bs[i,j] = -D[i,m+j]
        end
        bs[m+j,j] = one(T)
    end
end

# ガウス消去法
function _gauss_elimination!(D::AbstractMatrix{T}) where T
    m, n = size(D)

    # 一列目の行基本変形
    D[1,:] ./= D[1,1]
    for i in 2:m
        c = -D[i,1]
        for j in axes(D,2)
            D[i,j] += c * D[1,j]
        end
    end

    # 二列目以降の行基本変形
    i = j = 1
    while i <= m && j <= n
        if D[i,j] ≈ one(T)
            i += 1; j += 1; continue
        elseif D[i,j] ≈ zero(T)
            j += 1; continue
        else
            D[i,j:end] ./= D[i,j]
            for k in 1:i-1
                c = D[k,j]
                for l in i:n
                    D[k,l] -= c * D[i,l]
                end
            end
            i += 1
        end
    end
end
# シュミットの直交化
function _orthonormalize_by_schmidt!(bs::AbstractMatrix)
    for i in axes(bs,2)
        bi = @view bs[:,i]
        for j in 1:i-1
            bj = @view bs[:,j]
            corr = bj'bi # 複素数の場合内積の順序を間違えてはいけない
            bi .-= corr * bj
        end
        bi ./= sqrt(sum(abs2(x) for x in bi))
    end
end

# make dimension reduction matrix T
function _make_reduce_dim_mat(B, F, CIR, nfft, ntx, nrx)
    nrem = nrx - ntx
    L = size(CIR,1)-1 # チャネル長
    N = nfft # FFT数
    T = zeros(eltype(CIR), N*nrem, L*ntx) # 次元削減行列

    for k in 1:nrem
        for n in 1:ntx
            for l in 1:L
                f̂ = @view T[(k-1)*nfft+1:k*nfft,(n-1)*L+l]
                for m in 1:nrx
                    Bm = @views Diagonal(diag(B[(m-1)*nfft+1:m*nfft, (k-1)*nfft+1:k*nfft]))
                    for j in l:L
                        @views f̂ .+= CIR[j+1,m,n]*(Bm'*F[:,j-l+1])
                    end
                end
            end
        end
    end
    T
end


# SIMO detection(MRC)
function _simo_detection(z, U, D, CFR)
    α = _gain(CFR)
    return (D'*U' * z) ./ α
end

# 利得計算
function _gain(CFR)
    nfft, nrx, ntx = size(CFR)
    α = zeros(Float64, nfft, ntx)
    for i in eachindex(α)
        α[i] = sum(abs2(CFR[i,j]) for j in axes(CFR,2))
    end
    α
end

# MIMO検出
function _mimo_detection(y, H, dinds, method, N0=nothing, modt=nothing)
    nfft, nrx, ntx = size(H)
    nts = size(y,3) # time slots
    x = Array{eltype(y)}(undef, ntx, nfft, nts)

    if method == :ZF
        for k in axes(y,3) # time
            for j in dinds[k]
                x[:,j,k] = @views _zf(y[:,j,k], H[j,:,:])
            end
        end
    elseif method == :MMSE
        for k in axes(y,3) # time
            for j in dinds[k]
                x[:,j,k] = @views _mmse(y[:,j,k], H[j,:,:], N0)
            end
        end
    elseif method == :SIC || method == :VBLAST
        slicer = x -> slice(modt, x)
        for k in axes(y,3) # time
            for j in dinds[k]
                x[:,j,k] = @views _vblast(y[:,j,k], H[j,:,:], N0, slicer)
            end
        end
    elseif method == :MLD
        refs = getrefs(modt)
        for k in axes(y,3) # time
            for j in dinds[k]
                x[:,j,k] = @views _mld(y[:,j,k], H[j,:,:], refs)
            end
        end
    end
    return x
end

# MIMO-ZF検出
_zf(y, H) = inv(H'*H) * y

# MIMO-MMSE検出
_mmse(y, H, N0) = inv(H'*H + N0*I)*H'*y

# OSIC検出
function _vblast(y, H, N0, slicer)
    nrx, ntx = size(H)
    x = Array{eltype(y)}(undef, ntx)
    T = eltype(y)
    r = copy(y)
    inds = [1:ntx;]
    for st in 1:ntx
        if st == 1
            W = H*inv(H'H + N0*I)
            i = argmin([norm(x) for x in eachcol(W)])
            tmp = @views dot(W[:,i],r)
            deleteat!(inds, i)
        else
            i = inds[1]
            H_ = @views H[:,inds]
            w  = inv(H_*H_' + N0*I) * H_
            tmp = dot(w, r)
        end
        x[i] = slicer(tmp)
        r .-= H[:,i] * x[i]
    end
    #=
    H = H'*H

    w = Array{T}(undef,ntx)
    for (stage,i) = enumerate(ord)
        H_ = @view H[:,ord[stage:end]]
        hi = @view H[:,i]
        w .= inv(H_*H_' + N0*I) * hi
        x[i] = slicer(dot(w, r))
        if stage < ntx
            r .-= hi * x[i]
        end
    end
    =#
    x
end

_mld(y, H, refs) = MIMO.mld(y, H, refs)

# QR分解によるMIMOにおけるGSCの実装
function equalize_by_qr(ofdm::Ofdm, rx_frame::AbstractArray{T}, cir, cfr, N0; detection=:SIC, approx=false, option=()) where T
    nfft, ngi, nTs = ofdm.nfft, ofdm.ngi, ofdm.nTs
    ntx, nrx = ofdm.ndim
    nsym = size(rx_frame, 2)

    # 受信信号のstack
    r = Array{T}(undef, nfft*nrx, nsym)
    _stackvec!(r, rx_frame)

    # MMSEフィルタの生成
    W, R, order, N0 = _make_mmse_filter_by_qr(cfr, cir, N0, nfft, ngi, (ntx,nrx), detection, approx)

    # 信号検出
    y = Array{T}(undef, ntx, nfft, nsym)
    _destackvec!(y, W'*r)
    dinds = get_idxs_list(ofdm, :data)
    rxdata = _mimo_detection_by_qr(y, R, order, dinds, detection; option...)
    return @views rxdata[:, get_all_idxs(ofdm, :data)]
end

# 行列のstack
function _stackmat(H::AbstractArray{T,4}) where T
    N, Nr, Nt = size(H,1), size(H,3), size(H,4)
    Hs = permutedims(H, [3,1,4,2])
    Hs = reshape(Hs, Nr*N, Nt*N)
    Hs
end

# MMSEフィルタの生成
function _make_mmse_filter_by_qr(cfr, cir, N0, nfft, ngi, ndim, method, approx=false)
    # コヒーレント行列, ブロック行列の生成
    ntx, nrx = ndim
    # 干渉行列の生成・stack
    Hisi = Array{eltype(cfr)}(undef, nfft*nrx, nfft*ntx)
    Hisi_, _ = gen_coef_matrix(cir, nfft, ngi, ndim, domain=:freq)
    _stackmat!(Hisi, Hisi_)
    D, B, R, ord = _make_filter_matrix_by_qr(cfr)
    Rin = 2*Hisi*Hisi' + N0*I            # 干渉共分散行列

    # MMSEフィルタ行列生成
    if approx
        F = dftmtx(nfft, normalize=true)
        T = _make_reduce_dim_mat(B, F, cir, nfft, ntx, nrx)
        BT = B*T
        V  = I - (BT * inv(BT'*Rin*BT) * BT' * Rin)
    else
        U = inv(B'*Rin*B)*B'*Rin
        V = I - B*U
        #V = (I - B * inv(B'*Rin*B) * B'* Rin)
    end
    W = V * D
    if method == :MMSE || method == :OSIC
        N0 = tr(V'Rin) / (ndim[2]*nfft)
    end
    W, R, ord, N0
end



# フィルタ生成のための行列生成
function _make_filter_matrix_by_qr(CFR)
    nfft, nrx, ntx = size(CFR,1), size(CFR,2), size(CFR,3)
    rdim = nrx-ntx # redundant dimension

    # memory allocate
    D = spzeros(eltype(CFR), nrx*nfft, ntx*nfft)
    B = spzeros(eltype(CFR), nrx*nfft, nfft*rdim)
    R = zeros(eltype(CFR), ntx, ntx, nfft)
    orderlist = zeros(Int64, ntx, nfft) # column sort order

    # コヒーレント行列Dおよびブロック行列Bの生成
    for k = 1:nfft
        H = @view CFR[k,:,:]                             # チャネル行列
        orderlist[:,k] = columnsort(H, :SIR, rev=false)  # column sort
        H = H[:, orderlist[:,k]]
        Q_, R_  = qr(H)                                  # QR分解
        R[:,:,k] .= @views R_[1:ntx,1:ntx]               # 三角行列Rの保存

        for j = 1:ntx
            for i = 1:nrx
                D[k+nfft*(i-1), k+nfft*(j-1)] = Q_[i,j]
            end
        end
        for j = 1:rdim
            for i = 1:nrx
                B[k+nfft*(i-1), k+nfft*(j-1)] = Q_[i, end-rdim+j]
            end
        end
    end
    D, B, R, orderlist
end


# 信号検出
function _mimo_detection_by_qr(y, R, orderlist, dinds, method=:SIC; modt=nothing, M=0)
    x = similar(y)
    if method == :ZF
        for j in axes(y,3) # 時間軸
            for i in dinds[j]
                @views _back_subs!(x[orderlist[:,i],i,j], y[:,i,j], R[:,:,i])
            end
        end
    elseif method == :SIC
        refs = getrefs(modt)
        slicer = x -> slice(modt, x)
        for j in axes(y,3) # 時間軸
            for i in dinds[j]
                @views _sic!(x[orderlist[:,i],i,j], y[:,i,j], R[:,:,i], slicer)
            end
        end
    elseif method == :QRMLD
        M ==0 && error("パス数Mを指定してください。")
        refs = getrefs(modt)
        for j in axes(y,3)
            for i in dinds[j]
                @views _qrmld!(x[:,i,j], y[:,i,j], R[:,:,i], orderlist[:,i], refs, M)
            end
        end
    elseif method == :OSIC
        slicer = x -> slice(modt, x)
        for k in axes(y,3) # time
            for j in dinds[k]
                x[:,j,k] = @views _osic(y[:,j,k], H[j,:,:], N0, slicer)
            end
        end
    else
        error()
    end
    x
end


# ZF: 後退代入
function _back_subs!(x, y, R)
    for i in length(y):-1:1
        y_ = y[i]
        for j in i+1:length(y)
            y_ -= R[i,j]*x[j]
        end
        x[i] = y_ / R[i,i]
    end
end

# SIC
function _sic!(x, y, R, slicer)
    ntx = length(x); metric = 0.0;
    for i in ntx:-1:1
        tmp = y[i]
        for j in i+1:ntx
            tmp -= R[i,j] * x[j]
        end
        tmp /= R[i,i]
        x[i] = slicer(tmp)
        #=
        ri = y[i]; metric = 1e5; ind = 0
        # 干渉除去
        i < ntx && for j in i+1:ntx
            ri -= R[i,j]*x[j]
        end
        for j in eachindex(refs)
            tmp = abs2(ri - R[i,i]*refs[j])
            if metric>tmp
                metric = tmp
                ind = j
            end
        end
        x[i] = refs[ind]
        =#
    end
end

# QRMLD
function _qrmld!(x, z, R, ord, refs, M)
    x .= MIMO.qrmld(z, R, ord, refs, M)[1]
    return
end

end # module
