# Turbo.jl
module Turbo

using Random
using ..Conv: Trellis, APPDecoder
using ..FEC: AbstractEncoder, AbstractDecoder

export TurboEncoder,
       TurboDecoder,
       turboenc,
       turbodec


# Turboエンコーダ
struct Encoder <: AbstractEncoder
    enc1::Trellis            # トレリス1
    enc2::Trellis            # トレリス2
    rate::Rational{Int64}    # 符号化率
    bitlen::Int              # 入力ビット長
    tailbitlen::Int64        # 終端ビット数
    permind::Vector{Int64}   # 並び替え
    isterm::Bool             # トレリス終端の有無
    ispunc::Bool             # パンクチチャリング(true/false)
    P::Union{Nothing,Array{Bool}} # パンクチャリング行列
end
function Encoder(poly; fbpoly=poly[1], bitlen, isterm=true, rate=1//3)
    trellis1  = poly2trellis(poly[1:1, 1:2]; fbpoly=fbpoly, isterm=isterm)
    trellis2  = poly2trellis(poly[1:1, 1:2:3]; fbpoly=fbpoly, isterm=isterm)
    permind   = randperm(MersenneTwister(1234), bitlen)
    n_tailbit = 2*(maximum(trellis1.K) - 1) + 2*(maximum(trellis2.K) - 1)
    P, ispunc = rate==1//3 ? (nothing,false) : (get_punc_mtx(inout, rate),true)
    Encoder(trellis1, trellis2, rate, bitlen, n_tailbit, permind, isterm, ispunc, P)
end

# Turbo符号化(functor)
function (e::Encoder)(input)
    # RSC1, RSC2符号化
    rsc_out1 = convenc(e.enc1, input, isterm=d.isterm)
    rsc_out2 = convenc(e.enc2, view(input,e.permind) , isterm=e.isterm)
    parity1 = @view rsc_out1[2:2:length(input)*2]
    parity2 = @view rsc_out2[2:2:length(input)*2]

    # RSC1+RSC2 -> パンクチャ
    enc_bits = vcat(input', parity1', parity2') |> vec
    e.ispunc && (enc_bits=puncture(enc_bits, d.P))

    # 終端用ビット付加
    if e.isterm
        tail1 = view(rsc_out1, 2*e.n_bits+1:length(rsc_out1))
        tail2 = view(rsc_out2, 2*e.n_bits+1:length(rsc_out2))
        tailbits = vcat(tail1, tail2)
        enc_bits = append!(enc_bits, tailbits) # tailbit付加
    end
    enc_bits
end



# Turboデコーダ
struct Decoder <: AbstractDecoder
    dec1::APPDecoder             # APPdecoder1
    dec2::APPDecoder             # APPdecoder2
    rate::Rational{Int64}        # 符号化レート
    bitlen::Int64                # 入力ビット長
    tailbitlen::Int64            # 終端ビット数
    permind::Vector{Int64}       # 並び替えindex
    isterm::Bool                 # トレリス終端
    ispunc::Bool                 # パンクチチャ
    niter::Int64                 # 繰り返し回数
    method::String               # 復号方式
end

function Decoder(enc::Encoder; niter=6, method="max")
    decoder1 = APPDecoder(enc.enc1; n_bits=enc.n_bits, method=method, isterm=enc.isterm) # 事後確率復号器1
    decoder2 = APPDecoder(enc.enc2; n_bits=enc.n_bits, method=method, isterm=enc.isterm) # 事後確率復号器2
    Decoder(decoder1, decoder2, enc.rate, enc.n_bits, enc.n_tailbit, enc.permind, enc.isterm, enc.ispunc, enc.P, niter, method)
end

"""
Turbo復号メソッド
"""
function (d::Decoder)(Lch; target=:Input)
    T = eltype(Lch)
    ntb = div(d.n_tailbit, 4) # 各要素符号化器の終端ビット数
    bitlen = d.n_bits + ntb # 情報ビット数+終端ビット数
    if d.isterm
        tailbits = Lch[end-d.tailtbitlen+1:end] # tail-bits
        tailbits = reshape(tailbits, 2, :, 2) # tailbits
        Lch = Lch[1:length(Lch)-d.tailbitlen] # 情報ビット＋検査ビット
    end
    if d.ispunc
        output = zeros(eltype(Lch), 3, d.n_bits) # output array
        Lch = depuncture(output, Lch, d.P) # depuncture
    else
        Lch = reshape(Lch, 3, :)
    end

    # 配列確保
    if target == :Input
        Lapp = zeros(T, 1, bitlen) # 事後値値配列
    else
        Lapp = zeros(T, 3, bitlen) #
    end
    Lext = zeros(T, 1, d.n_bits)  # 事前値配列
    @views Lch1 = hcat(Lch[1:2, :], tailbits[:,:,1])
    @views Lch2 = hcat(vcat(Lch[1:1, d.permind], Lch[3:3,:]), tailbits[:,:,2])

    # 繰り返し復号
    for iter in 1:d.niter
        if iter < d.niter || target==:Input
            La1 = Lext
            bcjrdec!(Lapp, Lch1, La1; params=d.dec1, target=:Input)
            @views _subtract(Lext[1,:], Lapp[1,1:bitlen], La1[1,1:bitlen])
            La2 = Lext
            bcjrdec!(Lapp, Lch2, La2; params=d.dec2, target=:Input)
            @views _subtract(Lext[1,:], Lapp[1,1:bitlen], La2[1,1:bitlen])
        else
            @views bcjrdec!(Lapp[1:2,:], Lch1, La1; params=d.dec1, La=La1, target=:Output)
            @views _subctract(Lext, Lapp, La1)
            La2 = Lext
            @views bcjrdec!(Lapp[[1,3],:], Lch2, La2; params=d.dec2, La=La2, target=:Output)
        end
    end

    # 出力
    if target==:Input
        Lout = Lapp
        if d.isterm
            Lout = Lout[1,1:d.n_bits,2] # 情報ビット配列
            Lout[d.permind] = Lout # デインターリーブ
        end
        return Lout

    elseif target==:Output
        tailbits = vec(Lapp[:,end-ntb+1:end,:])
        # 出力配列
        @views Lapp = vcat(Lapp[1:1,1:n_bits,2], Lapp[2:2,1:n_bits,1], Lapp[2:2,1:n_bits,2])
        # デインタリーブ
        Lapp[1,d.permind] = Lapp[1,1:n_bits]
        Lout = vec(Lapp)
        append!(Lapp, tailbits)
        if d.puncturing # パンクチャリング
            Lout = puncture(Lout, d.P) # puncture
        end
        return Lout
    end
end

function _subtract!(Lext, Lapp, La)
    for n in eachindex(Lext)
        Lext[n] = Lapp[n] - La[n]
    end
end


# alias
TurboEncoder = Encoder
TurboDecoder = Decoder



end # module
