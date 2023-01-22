# kalman_equalization.jl

function kalman_equalization(ofdm::Ofdm, r::AbstractArray{T}, h, N0, inds, slice; n_iter) where T
    r = reshape(r, ofdm.n_Ts, :, ofdm.ndims[2])
    N, G, M = ofdm.nfft, ofdm.ngi, size(r,2)
    L = size(h,1)-1
    y = similar(r)                               # 等化受信信号(時間領域)
    Y = Array{T}(undef, N, M)                    # 等化受信信号(周波数領域)
    X = zeros(T, N, M)                           # 判定シンボルの格納配列
    x = Array{T}(undef, N+G, M)
    # x = reshape(x, N+G, :)
    FFT  = @views plan_fft(x[G+1:end,1], 1)
    IFFT = @views plan_ifft(x[G+1:end,1], 1)

    # チャネル推定
    H  =  OFDM.ChEst.convert_to_cfr(h, N, G)
    h = reshape(h, :)

    # カルマンフィルタパラメータセットアップ
    n_state_dim = length(h)
    n_obs_dim   = ofdm.ndims[2]
    params = _setup(T, h, N0, (N+G)*M, n_state_dim, n_obs_dim)

    _kalman_equalization!(x, X, y, Y, r, h, H, FFT, IFFT, N, G, L, params, inds, slice, n_iter)
    X
end

# セットアップ
function _setup(T, h, N0, len, n_state_dim, n_obs_dim)
    x = zeros(T, n_state_dim, len)                  # 状態
    V = zeros(T, n_state_dim, n_state_dim, len)     # 状態の共分散行列
    G = zeros(T, n_state_dim)
    H = reshape(h, 1, :)
    Q = zeros(T, n_state_dim, n_state_dim)
    R = zeros(T, n_obs_dim, n_obs_dim)

    # 初期化
    Q = 1.0         # システム雑音分散の初期値
    R = N0            # 雑音分散
    x, V, H, R, Q
end


function _kalman_equalization!(x, X, y, Y, r, h, H, FFT, IFFT, N, G, L, params, inds, slice, n_iter)
    iter = 0
    while iter < n_iter
        y .= r
        if iter == 0
            for k in 1:size(r,2)
                # ISI subtraction
                if 1<k
                    x[G+1:end,k-1] .= @views (IFFT * X[:,k-1]) * sqrt(N)
                    @views _isi_remove!(y[:,k], x[:,k-1], h, G, L)
                end

                # Kalman filtering & ISI,ICI compensation
                @views _kalman_filter!(y[:], params...)
                Base.print_matrix(stdout, params[2][1,1,:])
                x[G+1:end,k] .= @view params[1][1,G+1:end]

                # ICI compensation
                @views _ici_remove!(y[:,k], x[:,k], h, G, L)

                # Data detection
                Y[:,k] = (FFT * @view y[G+1:end,k]) / sqrt(N)
                @views _data_detection!(X[:,k], Y[:,k], H, inds[:,k], slice)
            end
        else
        end
        iter += 1
    end
end



# カルマンフィルタ
function _kalman_filter!(y, x, V, H, R, Q)
    x_ = x[:,1]
    V_ = V[:,:,1]
    x_[:]   .= zero(eltype(x))
    V_[:,:] .= zero(eltype(V))
    V_[diagind(V_)] .= 10000.0

    for n in 17:size(y,1)
        @views _filter!(y[n], x[:,n], V[:,:,n], x_, V_, H, R)
        @views _predict!(x_, V_, x[:,n], V[:,:,n], Q)
    end
end


# フィルタリングステップ
function _filter!(y, x, V, x_, V_, H, R)
    K = (V_[1] * H'[1]) * inv(real(H*V_*H')[1] + R[1])
    # x .= x_ + K * (y - h' * x_)
    x .= x_; V .= V_
    x[1] = K * (y - (H*x_)[1])
    V[1] = (1.0 - K*H[1])*V_[1]
    # V[diagind(V)] .= real(diag(V))
end

# 予測ステップ
function _predict!(x_, V_, x, V, Q)
    # next state
    x_[1] = 0.0
    for i in 2:length(x)
        x_[i] = x[i-1]
    end

    # next covariance
    m = size(V,1)
    V_[2:m,2:m] .= @view V[1:m-1,1:m-1]
    V_[1] = Q[1]
end
