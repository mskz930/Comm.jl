# mld.jl
using LinearAlgebra: dot, norm, Diagonal, diag, diagm, diagind, I
using InvertedIndices: Not
using Comm.DigitalModulation

# methods
include("sorts.jl")

struct MLD{T,U}
  refs::Vector{T}
end
function MLD(modt)
  refs = get_refs(modt)
  MLD{eltype(refs)}(refs)
end

function mld!(x, y, H, refs)
  nrx, ntx = size(H)
  x, ml_list, ml_metric = _mld!(x, y, H, refs)
  x, ml_list, ml_metric
end

function mld(y::AbstractVector{T}, H::AbstractMatrix{T}, refs) where {T}
  x = zeros(T, size(H, 2))
  x, list, metric = mld!(x, y, H, refs)
  return x, list, metric
end

# 最尤パス探索
function _mld!(x, y, H, refs)
  # ローカル変数
  nref = length(refs)
  x_tmp = similar(x) # candidate vector
  list_tmp = ones(Int32, length(x_tmp))
  ml_metric = Inf
  ml_list = copy(list_tmp)

  M = length(refs)
  n = length(x_tmp)
  # MLメトリックの計算
  for k = 1:M^n
    k > 1 && _next_list!(list_tmp, nref)
    _candidate!(x_tmp, list_tmp, refs) # 候補ベクトルを取得
    metric_tmp = @views _calc_metric(y, H, x_tmp) # メトリック計算
    # @show n, smetric_tmp
    if ml_metric > metric_tmp
      ml_list .= list_tmp
      ml_metric = metric_tmp
    end
  end
  # @show ml_metric
  x_ml = @views refs[ml_list]
  x_ml, ml_list, ml_metric
end

# インデックスから候補ベクトルを取得
function _candidate!(x_tmp, tmp_list, refs)
  for i in eachindex(x_tmp)
    x_tmp[i] = refs[tmp_list[i]]
  end
end

# 次リスト(index)の計算
function _next_list!(list, modN)
  i = 0
  flag = true
  while flag
    list[end-i] = list[end-i] % modN + 1
    flag = list[end-i] == 1
    i += 1
  end
end

# リストの保存
function _save_list!(list, tmp_list, j)
  for i in axes(list, 1)
    list[i, j] = tmp_list[i]
  end
  list
end

# MLメトリックの計算
function _calc_metric(y, H, x_tmp)
  metric = 0.0
  @inbounds for i in axes(H, 1)
    tmp = zero(eltype(y))
    for j in axes(H, 2)
      tmp += H[i, j] * x_tmp[j]
    end
    metric += abs2(y[i] - tmp)
  end
  metric
end
