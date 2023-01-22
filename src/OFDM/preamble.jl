# preamble.jl

"""
Preambleパラメタ
"""
struct Preamble
    type::String
    data::Array{ComplexF64}
    size::Tuple{Int64,Int64}
    Preamble(type, data) = new(type, data, size(data))
end

# プリアンブルシンボルの生成
function gen_preamble(info::Preamble, n_fft, n_gi; kind=:long)
    pilot_symbols = randqam(n_fft, M=4)
    preamble = ifft(pilot_symbols)
    preamble = [preamble[end-n_gi+1:end], preamble] # CP付加
end
