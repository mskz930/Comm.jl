
# ビタビ復号(viterbi decoding)
function vitdec(r, trellis::Trellis, dectype="hard")
    # 配列確保
    x = BitArray(undef, bitlen)
    metrics = Array{Float64}(undef, trellis.n_state, bitlen)
    nextstates = trellis.nextstates
    if dectype=="hard"
        viterbi_hard(x, y, nextstates, outputs)
    end
    viterbi!(x, y, nextstates, outputs)
end

# 硬判定ビタビ復号
function viterbi!(x, y, numstates, nextstates)
    # ビタビ復号
    for n in 2:size(y,2)
        # メトリック初期化
        @views _init(metric[:,n])

        # パスメトリック計算

        # メトリックの更新

    end
end

# メモリ初期化
function _init(metric)
    MINF = 100000
    @inbounds @simd for i in eachindex(memory)
        metric[i] = MINF
    end
end


# トレースバック
function traceback()
end
