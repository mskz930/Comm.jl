

# フレーム情報
struct DataFrame{T,N}
    A::Array{T,N}
    idx_table::Array{Int8,N}
end
function DataFrame(T::DataType, n_data::Int, n_fft, n_Tx, pilot, idxs, max_frame_len=Inf)
    frame_len= getframesize(n_data, n_Tx, pilot, idxs) # フレーム長
    A   = zeros(ComplexF64, n_fft, frame_len, n_Tx)
    idx_table = zeros(Int8, size(A)...)
    DataFrame{T,3}(A, idx_table)
end


Base.show(io::IO, f::DataFrame) = print(io, "$(typeof(f))")
Base.show(io::IO, m::MIME"text/plain", f::DataFrame) = print(io, """$(typeof(f)):
    idxs => $(typeof(f.idxs)): with size $(size(f.idxs))
    idxs_lis =>
    on   => $(f.on)
""")

# DataFrame.idxsを0初期化する
function _init_DataFrame_idxs!(DataFrame::DataFrame)
    DataFrame.idxs .= 0
    DataFrame.idxs_list = Dict{Symbol,Vector{Vector{Int8}}}(
        :data  => Vector{Int8}[],
        :pilot => Vector{Int8}[],
        :guard => Vector{Int8}[],
    )
    return
end

# 割り当て済みか確認
isallocated(f::DataFrame) = return (f.on[:data] + f.on[:pilot]) > 0

"""
    init_DataFrame!(DataFrame1, DataFrame2)
フレームを初期化する
    DataFrame1 <= DataFrame2
"""
function init_DataFrame!(DataFrame::DataFrame, DataFrame1::AbstractArray{T}, DataFrame2::AbstractArray{U}) where {T<:SoftSymbol, U}
    for j = axes(DataFrame2,2)
        for i = axes(DataFrame2,1)
            if DataFrame.idxs[i,j] != -1
                for k = axes(DataFrame2,3)
                    DataFrame1[i,j,k] = SoftSymbol(DataFrame2[i,j,k], 0.0)
                end
            end
        end
    end
    DataFrame1
end

"""
    get_idxs(f::DataFrame, target)

OFDMフレームのインデックスを取得
    DataFrame  : DataFrame object
    target : {:data, :pilot, :dummy, :all}
"""
function get_idxs(f::DataFrame, target)
    if target == :data
        return findall(x -> x == -1, f.idxs)
    elseif target==:pilot
        return findall(x -> x > 0, f.idxs)
    elseif target == :dummy
        return findall(x -> x == -2, f.idxs)
    elseif target == :all
        return findall(x -> x < 0, f.idxs)
    end
end

"""
    get_idxs_list(f::DataFrame, target, rng)

    DataFrame  : DataFrame Object
    target :
    rng    : range
"""
function get_idxs_list(f::DataFrame, target::Symbol, rng=axes(f.idxs,2))
    if target == :data
        @views [findall(x -> x==-1, f.idxs[:,j]) for j in rng]
    elseif target == :all_data
        @views [findall(x -> x<0, f.idxs[:,j]) for j in rng]
    elseif target == :pilot
        @views [findall(x -> x>0, f.idxs[:,j]) for j in rng]
    elseif target == :guard || target == :null
        @views [findall(x -> x==0, f.idxs[:,j]) for j in rng]
    elseif target == :all
        @views [LinearIndices(f.idxs[:,j]) for j in rng]
    end
end
function get_idxs_list(f::DataFrame, target::Integer, rng=axes(f,2))
    @assert rng[end] <= size(f,2)
    if target == -1
        return [findall(x->x<0, f.idxs[:,j]) for j in rng]
    elseif target > 0
        [findall(x->x==id, f.idxs[:,j]) for j in rng]
    elseif target == 0
        [findall(x->x==0, f.idxs[:,j]) for j in rng]
    else
        Int[]
    end
end
