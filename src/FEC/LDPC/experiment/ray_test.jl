# decoding_test.jl AWGN channel BER test
module Test

using JLD2, FileIO, Distributions
include("../LdpcCode_kai.jl")

src = string(pwd(), "/Mymodule/Commun/ErrorCorrecting/ldpc/code_set/PEG_irr10_ray_8192_16384.jld2")
G, H = load(src, "G" ,"H") # 生成行列, 検査行列読み込み
rate = size(G,2)//size(G,1) # 符号化レート
nbits = size(G,2)
ldpcenc = LDPC.LdpcEncoder(G)
ldpcdec = LDPC.LdpcDecoder(H)

function main(ebnodBs)
    ebnodBs = Vector(ebnodBs)
    ber_results = zeros(size(ebnodBs))
    for i in eachindex(ebnodBs)
        ber_results[i] = bercalc(ebnodBs[i])
        @show ebnodBs[i], ber_results[i]
        if ber_results[i] == 0
            break
        end
    end
    ber_results
end

function bercalc(ebnodB)
    m_ebnodB = ebnodB + 10*log10(rate) # EbN0
    N0 = 10^(-m_ebnodB/10) # SNR
    numbits=0; errors=0
    @time while (numbits<4e6) && (errors<1e4)
        tx_bits = rand(nbits) .> 0.5 # 送信ビット列
        enc_bits = ldpcenc(tx_bits) # 送信符号語
        tx_symbols = 2*enc_bits .- 1 # BPSKシンボル
        h = abs.(randn(ComplexF64, length(tx_symbols))) # チャネル係数
        # h = abs.(randn(ComplexF64, 1, Int(length(tx_symbols)/64)))
        # h = repeat(h, 64, 1)[:]
        rx_symbols = h .* tx_symbols #
        rx_symbols .+= sqrt(N0/2)*randn(Float64, size(tx_symbols)) # 受信シンボル
        Lch = 4.0*(rx_symbols.*h)./N0 # チャネル尤度
        dec_bits = ldpcdec(-Lch, 80; outtype=:hard) # sumproduct復号
        rx_bits = view(dec_bits,1:nbits)
        errors += sum(tx_bits .!== rx_bits) # 誤りビット数
        error = sum(tx_bits .!== rx_bits)
        numbits += nbits
    end
    ber = errors/numbits
    return ber
end
function test(ebnodB)
    m_ebnodB = ebnodB + 10log10(rate) # EbN0
    N0 = 10^(-m_ebnodB/10) # SNR
    tx_bits = rand(nbits) .> 0.5 # 送信ビット列
    enc_bits = ldpcenc(tx_bits) # 送信符号語
    tx_symbols = 2*enc_bits .- 1 # BPSKシンボル
    h = abs.(randn(Float64, 1, Int(length(tx_symbols)/64))) # チャネル係数
    h = repeat(h,64,1)
    rx_symbols = h .* tx_symbols #
    rx_symbols .+= sqrt(N0/2)*randn(Float64, size(tx_symbols)) # 受信シンボル
    Lch = 4.0*rx_symbols./N0 # チャネル尤度
    @time dec_bits = ldpcdec(-Lch, 50; outtype=:soft) # sumproduct復号
    dec_bits = -dec_bits .> 0
    rx_bits = view(dec_bits,1:nbits)
    sum(tx_bits .!== rx_bits) # 誤りビット数
end

end
