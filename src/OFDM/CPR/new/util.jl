
"""
干渉成分行列Hisi, Hiciを計算する
"""


"""
SINR計算を計算
  - 各受信アンテナにおける平均SINRを計算する
"""
function sinr_calc(cir, N0, nfft)
    L, n_rx, n_tx = size(cir)
    sinrs = zeros(Float64,length(range)) # 配列
    n_path = L - 1
    Ps = Pisi = Pici = 0.0
    for q = 1:n_rx # n_rx
      for p = 1:n_tx # n_tx
        for i in 1:L
            
        end
    end
    return sinrs
end
