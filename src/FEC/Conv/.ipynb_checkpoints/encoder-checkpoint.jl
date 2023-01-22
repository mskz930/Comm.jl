
"""
畳み込み符号化(Encoding)
    convenc(input, trellis)
    convenc!(output, input, trellis)
    (::Trellis)(input; kwargs...)
    (::Trellis)(output, input; kwargs...)

"""
function encode(input; trellis::Trellis=trellis, kwargs...)
    puncture = get(kwargs, :puncture, false) # パンクチャの処理

    # 前処理
    N = length(input) # 入力ビット長
    L = maximum(trellis.K) # 拘束長
    M = N + isterm(trellis)*(L - 1) # 出力長
    if ndims(input)==1
        input = reshape(input, trellis.size[1], :)
    end

    # 畳み込み符号化
    output = zeros(Bool, trellis.size[2], M)
    output = encode!(output, input; trellis=trellis, puncture)
    return vec(output)
end
function encode!(output, input; trellis::Trellis=trellis,  puncture=NamedTuple())
    K, N = size(input)  # 入力次元, 入力長
    state = 0 # 初期状態
    base = 0
    n = 1
    while n <= N
        is = @views bin2dec(input[:,n])       # 入力シンボル
        os = trellis.outputs[state+1, is+1]   # 出力シンボル
        @views dec2bin!(output[:,n], os)
        state = trellis.nextstates[state+1, is+1] # 次状態
        n += 1
    end
    while trellis.terminate && state > 0 && n <= size(output,2)
        is = @views argmin(trellis.nextstates[state+1, :]) - 1
        os = trellis.outputs[state+1, is+1]
        @views dec2bin!(output[:,n], os)
        state = trellis.nextstates[state+1, is+1]
        n += 1
    end
    trellis.terminate && state!=0 && error("正常に終端できませんでした．state=$state")
    output
end

convenc = encode
convenc! = encode!

(t::Trellis)(input; kwargs...) = encode(input, t; kwargs...)
(t::Trellis)(output, input; kwargs...) = encode!(output, input, t; kwargs...)


struct Encoder
end

function (enc::Encoder)(x)
end
