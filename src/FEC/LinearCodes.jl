# 線形符号用モジュール
module LinearCodes

using Binary
using LinearAlgebra: rank
using SparseArrays: SparseMatrixCSC, AbstractSparseMatrix
import LinearAlgebra: rank

include("crc.jl")

export Generator, ParityCheck

struct Generator{T <: AbstractMatrix}
    shape::Tuple{Int,Int}
    P::T

    function Generator(G::AbstractMatrix)
        n,k = shape = size(G)
        P = G[k+1:end, :]
        new{typeof(G)}(shape, P)
    end
end

rank(G::Generator) = G.shape[1] - G.shape[2]



struct ParityCheckMatrix{T<:AbstractMatrix}
    mat::T
    nodelist::Dict{:Symbol, Vector{Int}}

    function ParityCheck(H::T, rank) where T <: AbstractMatrix
        @assert maximum(H) <= 1
        if typeof(H) <: DenseArray
            H = sparse(H)
        end
        new{T}(H, rank)
    end
end

# 整数のバイナリに対する排他的論理和
function xorsum(s::T) where T
    r = zero(T)
    while s>0
        r ⊻= s & 1
        s >>= 1
    end
    r
end

# num of bit of "one"
function bin_sum(d::T) where T
    sum = zero(T)
    while d>0
        sum += d & 1
        d >>= 1
    end
    sum
end

"""
ハミング距離
"""
hamming(x::T, y::T) where T <: Integer = bin_sum(x ⊻ y)

"""
ユークリッド距離
"""
euclidean(x::T, y::T) where T <: AbstractFloat = abs(x .- y)

""" get_matrix: 生成行列の取得 """
function get_genmatrix(H::AbstractArray)
    rref_H, pivot, rank = rref(H)
    pivot_else = setdiff(1:size(H,2), pivot)
    k = length(pivot_else) # 情報ビット数
    P = view(rref_H, :, pivot_else) # 線形従属な部分
    eye = Matrix(I,k,k) # 単位行列
    G = [eye; P] # 生成行列
    H = H[:,vcat(pivot_else,pivot)] # パリティ検査行列
    return G, H
end


""" rref!: 行基本変形(reduced-row-echelon-form) """
function rref!(A::AbstractMatrix{Bool})
    M, N = size(A) # 行, 列
    row = 1; col = 1 # 初期値
    pivot = Int64[]
    while(row<=M) && (col<=N)
        # col列のrow以下の非ゼロ行を取得
        # index = findall(x->x==true,A[:,col])
        a = view(A,:,col)
        index = LinearIndices(a)[a]
        lower_index = filter(x -> x>row, index) # 現在の行よりも下の行インデックス
        if (A[row,col] == false) && (!isempty(lower_index))
            # 行の入れかえが可能であるとき
            i = lower_index[1] # 最も近い行インデックス
            index = filter(x-> x!==i, index)
            A_ = A[row,:]
            A[row,:] = view(A,i,:) # 行入れ替え
            A[i,:] = A_
        elseif (A[row,col] == false) && (isempty(lower_index))
            # 行の入れかえが不可能なとき
            col += 1
            continue # 右の列に移動する
        end
        # col列を単位行列にする
        push!(pivot,col)
        b = view(A,row,:)
        colinds = LinearIndices(b)[b]
        for j in colinds
            for i in index
                if i !== row
                    A[i,j] ⊻= A[row,j] # 各行にXOR加算
                end
            end
        end
        col += 1 # 右列に移動
        row += 1 # 下の行に移動
    end
    rank = length(pivot) # 行列のランク
    return A, pivot, rank
end
rref(A) = rref!(copy(A))

function rref!(A::AbstractMatrix{<:Float64})
    M, N = size(A) # 行, 列
    row = 1; col = 1 # 初期値
    pivot = Int64[]
    while(row<=M) && (col<=N)
        # col列のrow以下の非ゼロ行を取得
        # index = findall(x->x==true,A[:,col])
        a = view(A,:,col)
        index = LinearIndices(a)[a .> 0]
        lower_index = filter(x -> x>row, index) # 現在の行よりも下の行インデックス
        if (A[row,col] == 0) && (!isempty(lower_index))
            # 行の入れかえが可能であるとき
            i = lower_index[1] # 最も近い行インデックス
            index = filter(x-> x!==i, index)
            A_ = A[row,:]
            A[row,:] = view(A,i,:) # 行入れ替え
            A[i,:] = A_
        elseif (A[row,col] == 0) && (isempty(lower_index))
            # 行の入れかえが不可能なとき
            col += 1
            continue # 右の列に移動する
        end
        # col列を単位行列にする
        push!(pivot,col)
        colinds = findall(x->x > 0, view(A,row,:))
        for j in colinds
            for i in index
                if i !== row
                    A[i,j] += A[row,j] # 各行にXOR加算
                    A[i,j] %= 2
                end
            end
        end
        col += 1 # 右列に移動
        row += 1 # 下の行に移動
    end
    rank = length(pivot) # 行列のランク
    return A, pivot, rank
end




end # module end
