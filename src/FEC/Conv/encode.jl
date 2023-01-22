

mutable struct ConvolutionalEncoder
    trellis::Trellis
    bit_length::Int
    code_length::Int
end
function ConvolutionalEncoder(trellis::Trellis; bit_length=nothing, code_length=nothing)
    if bit_length === nothing
        bit_length = code_length = -1
    end
    ConvolutionalEncoder(trellis, bit_length, code_length)
end
ConvEncoder = ConvolutionalEncoder

size(enc::ConvEncoder) = (enc.bit_length[], enc.code_length[])
set!(enc::ConvEncoder, property::Symbol, value::Any) = setproperty!(enc, property, value)

"""
畳み込み符号化(Encoding)
    convenc(input, trellis)
    convenc!(output, input, trellis)
    (::Trellis)(input; kwargs...)
    (::Trellis)(output, input; kwargs...)

"""
function encode(input; trellis::Trellis=trellis, ispunc=false)

    # 前処理
    N = length(input) # 入力ビット長
    L = maximum(trellis.K) # 拘束長
    M = N + isterm(trellis) * (L - 1) # 出力長
    if ndims(input) == 1
        input = reshape(input, trellis.size[1], :)
    end

    # 畳み込み符号化
    output = zeros(Bool, trellis.size[2], M)
    output = encode!(output, input; trellis=trellis, ispunc=ispunc)
    return vec(output)
end
# function convenc(encoder::ConvEncoder, input::BitVector)
#     @assert encoder.bit_length == length(input)
# end
function encode!(output, input; trellis::Trellis=trellis, ispunc=false)
    K, N = size(input)  # 入力次元, 入力長
    state = 0 # 初期状態
    base = 0
    n = 1
    while n <= N
        is = @views bin2dec(input[:, n]) # binary => decimal
        os = trellis.outputs[state+1, is+1] # next output
        @views dec2bin!(output[:, n], os) # decimal => binary
        state = trellis.nextstates[state+1, is+1] # next state
        n += 1
    end
    while trellis.terminate && state > 0 && n <= size(output, 2)
        is = @views argmin(trellis.nextstates[state+1, :]) - 1
        os = trellis.outputs[state+1, is+1]
        @views dec2bin!(output[:, n], os)
        state = trellis.nextstates[state+1, is+1]
        n += 1
    end
    trellis.terminate && state != 0 && error("正常に終端できませんでした．state=$state")
    output
end

convenc = encode
convenc! = encode!

(t::Trellis)(input; kwargs...) = encode(input, t; kwargs...)
(t::Trellis)(output, input; kwargs...) = encode!(output, input, t; kwargs...)


struct Encoder
    trellis::Trellis
    input_size::Vector{Int}
    Encoder(trellis::Trellis) = new(trellis)
    Encoder(trellis::Trellis, input_size::Int) = new(trellis, [input_size])
end

function get_input_size(enc::Encoder)
    if isdefined(enc, :input_size)
        enc.input_size[1] > 0 && return enc.input_size[1]
        error("input size must be positive value.")
    else
        error("input size is not defined.")
    end
end

function set_input_size(enc::Encoder, input_size::Int)
    if isdefined(enc, :input_size)
        enc.input_size[1] = input_size
    else
        enc.input_size = [input_size]
    end
end


