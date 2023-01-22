# APPDec_copy.jl テスト用ファイル
module APPDec

include("../../Utils/BinaryTools.jl")
using .BinaryTools
using Commun.FEC.Conv

export APPDecoder, Trellis, poly2trellis

# 対数の最小値定義
const MINF = -1e20


"""
    APPDecoder

APP復号器オブジェクト
"""
struct APPDecoder
    n_inputs::Int64    # 入力ビット数
    n_outputs::Int64   # 出力ビット数
    n_states::Int64    # 状態数
    n_bits::Integer    # ビット長
    transmat::Dict{Symbol, Array{Vector{Float64},2}} # 状態遷移配列(:input, :output)
    mask::Matrix{Bool}       # マスク
    alpha::Matrix{Float64}   # alpha
    beta::Matrix{Float64}    # beta
    gamma::Array{Float64,3}  # gamma
    isterm::Bool             # トレリス終端
    operator::String         # 作用素
end

# コンストラクタ
function APPDecoder(trellis::Trellis; n_bits, operator="max", isterm=true)
    @assert method in ["max", "max*", "logmap"]
    outputs = dec2bin.(trellis.outputs, len=2^(trellis.n_outputs-1)) # 出力シンボル(バイナリ)
    nextstates = trellis.nextstates # 次状態リスト
    n_states = trellis.n_states # 状態数
    n_tailbits = isterm ? maximum(trellis.K)-1 : 0 # 終端ビット数
    transmat, mask = _gen_transmat(nextstates, outputs, n_states) # 遷移行列の取得
    bitlen = n_bits + n_tailbits # 出力ビット長
    alpha = Array{Float64}(undef, n_states, bitlen+1) # alpha
    beta  = Array{Float64}(undef, n_states, bitlen+1) # beta
    gamma = Array{Float64}(undef, n_states, n_states, bitlen) # gamma
    APPDecoder(trellis.n_inputs, trellis.n_outputs, trellis.n_states, n_bits, transmat, mask, alpha, beta, gamma, isterm, method)
end


# 状態遷移行列生成
function _gen_transmat(nextstates::AbstractArray, outputs::AbstractArray, n_states::Integer, f = x -> 2x .- 1)
    outputs = map(f, outputs) # シンボル化
    insize = size(nextstates,2) # 入力サイズ
    input_mat = fill(Float64[], n_states, n_states) # 遷移行列(入力)
    output_mat = fill(Float64[], n_states, n_states) # 遷移行列(出力)
    transmat = Dict(:input=>input_mat, :output=>output_mat) # 遷移行列
    mask = zeros(Bool,n_states,n_states) # マスク
    for i in 1:n_states
        for j in 1:insize
            ns = nextstates[i,j] # 次状態
            # 出力ビット列を挿入
            transmat[:input][i,ns+1] = 2 .* dec2bin(j-1, len=Int(log2(insize))) .- 1
            transmat[:output][i,ns+1] = outputs[i,j]
            mask[i,ns+1] = 1
        end
    end
    transmat, mask
end



"""
    (dec)()
最大事後確率復号関数
"""
function (dec::APPDecoder)(Lin::T, La::Union{Nothing,T}=nothing; target="input", isterm=true, ispunc=false) where T<:AbstractArray
    # 入力データ整形
    if ispunc
        Lin = depunc(Lin, dec.P) # depuncture
    else
        Lin = reshape(Lin, dec.n_inputs, :) # reshape
    end
    # 出力配列確保
    bitlen = size(Lin, 2) # 情報ビット＋終端ビット
    n_rows = target==:input ? dec.n_inputs : dec.n_outputs
    Lout = zeros(Float64, n_rows, bitlen)
    # BCJR復号の実行
    _bcjrdec!(Lout, Lin, La, dec.alpha, dec.beta, dec.gamma, dec.transmat, dec.mask, target, isterm, max)
    # 出力データ整形
    if target==:input
        return Lout[1:dec.n_bits]
    else
        return vec(Lout)
    end
end


function (dec::APPDecoder)(Lout::T, Lin::T, La::Union{Nothing,T}=nothing; target="input", isterm=true, ispunc=false) where T<:AbstractArray
    # BCJR復号の実行
    _bcjrdec!(Lout, Lin; La=La, params=dec, target=target)
    # 出力データ整形
    if target==:input
        return Lout[1:dec.n_bits]
    else
        return vec(Lout)
    end
end



"""
    bcjrdec!()

BCJR復号
    arguments:
        Lin, La, Lout, insize: チャネル値, 事前値, 事後値
        insize, outsize, insymbol, outsymbol: 入力ビット数, 出力ビット数, 入力ビット配列, 出力ビット配列
        indices, n_states, outtype, logmap, isterm: 遷移行列のインデックス, 状態数, 出力尤度, 計算手法, 打ち切り情報
    returns:
"""
function _bcjrdec!(Lapp, Lch, La, alpha, beta, gamma, transmat, mask, target=:input, isterm=true, operator=max)
    # パラメータ設定
    colinds = [findall(x->x > 0, view(mask,i,:)) for i in 1:size(mask,1)] # 遷移行列の各行の非ゼロな列index
    inputs, outputs = transmat[:input], transmat[:output]
    # alpha, gamma, beta初期化
    _init!(alpha, beta, isterm)
    # 前向き計算
    _forward!(Lch, La, alpha, gamma, outputs, colinds, operator)
    # 後ろ向き計算
    symbol = target==:input ? insymbol : outsymbol
    _backward!(Lapp, alpha, beta, gamma, indices, symbol, operator)
end

# alpha, beta初期化
function _init!(alpha, beta, gamma, Lch, outputs, isterm)
    global MINF
    alpha[1,1] = 0.0; alpha[2:end,1] .= -MINF
    if isterm # トレリス終端の有無
        beta[1,end] = 0.0; beta[2:end,end] .= -MINF
    else
        beta[:,end] .= -log(size(alpha,1)) # 一様確率
    end
    # gamma 初期化
    for k in axes(gamma,3)
        for m in 1:n_outputs
            # 各出力シンボルに対する尤度を計算
            gamma_temp = 0.0
            for l in axes(Lch,1) # l-th bit channel LLR
                gamma_temp += Lch[l,k] * outputs[i,j,l]
            end
            # 同じ出力パスに対して代入
            for (i,j) in paths[m]
                gamma[i,j,k] = gamma_temp
            end
        end
    end
end

# 前向き計算
function _forward!(Lch, La, alpha, gamma, symbols, col_inds, operator)
    global MINF
    for n in axes(Lch,2)
        #= αの初期化
        for i in axes(alpha,1)
            alpha[i,n+1] = -MINF
        end
        =#
        alpha_max = -MINF
        for i in axes(gamma,1)
            # γの更新
            for j in col_inds[i]
                gamma_temp = 0.0
                for m in axes(Lch,1) # チャネル尤度加算
                    gamma_temp +=  (Lch[m,n]/2) * outputs[i,j][m]
                end
                if !isnothing(La) && n <= length(La)
                    for m in axes(La,1) # 事前尤度加算
                        gamma_temp += (La[m,n]/2) * symbols[i,j][m]
                    end
                end
                gamma[i,j,n] = gamma_temp
            end
            # α更新
            alpha_temp = -MINF
            for j in col_inds[i]
                alpha_temp = operator(alpha_temp, alpha[i,n]+gamma[i,j,n])
            end
            alpha_max  = alpha_max<alpha_temp ? alpha_temp : alpha_max # alphaの最大値を保存
        end
        # オーバーフロー対策(最大値を引く)
        for i in axes(alpha,1)
            alpha[i,n+1] -= alpha_max
        end
    end
end


# 後ろ向き計算
function _backward!(Lapp, alpha, beta, gamma, indices, symbol, func)
    global MINF
    for n in size(Lapp,2):-1:1
        # β初期化
        for i in axes(beta,1); beta[i,n] = -MINF; end

        # βの計算
        beta_max = -MINF
        for (i,j) in indices
            beta[i,n] = func(beta[i,n], gamma[i,j,n] + beta[j,n+1])
            beta_max = beta_max < beta[i,n] ? beta[i,n] : beta_max # alphaの最大値を探す
        end
        # オーバーフロー対策(最大値を引く)
        for i in axes(beta,1); beta[i,n] -= beta_max; end

        # 事後確率計算
        @views _appcalc!(Lapp[:,n], alpha[:,n], beta[:,n+1], gamma[:,:,n], indices, symbol, func, MINF)
    end
    Lapp
end

# 事後確率計算
function _appcalc!(Lapp, alpha, beta, gamma, indices, symbol, func, MINF)
    for k in axes(Lapp,1)
        A=-MINF; B=-MINF
        for (i,j) in indices
            if symbol[i,j][k] > 0 # if p(y|x=+1)
                A = func(A, alpha[i] + gamma[i,j] + beta[j])
            else # if p(y|x=-1)
                B = func(B, alpha[i] + gamma[i,j] + beta[j])
            end
        end
        # Lapp ~ max{p(x=+1|y)} - max{p(x=-1|y)}
        Lapp[k] = A-B
    end
end


end # module
