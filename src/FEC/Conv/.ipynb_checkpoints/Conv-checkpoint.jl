module Conv

using Binary
using SparseArrays
using UnPack
using DataStructures
using ..FEC: AbstractEncoder, AbstractDecoder

# source files
include("utils.jl")         # 汎用関数
include("trellis.jl")       # トレリス
include("encoder.jl")       # 符号化
include("vitdec.jl")        # ビタビ復号
include("decoder.jl")        # MAP復号

export convenc,   # 畳み込み符号化
       convenc!,
       puncture,  # パンクチャ
       depuncture, # デパンクチャ
       get_punc_mtx, # パンクチャ行列の取得
       Trellis, # トレリス型
       getrate,
       numofoutputs,
       make_transmat,
       poly2trellis, # トレリス型のコンストラクタ
       BCJRDecoder, # BCJR復号器パラメータ
       APPDecoder # 最大事後確率復号


end
