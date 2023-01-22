# 廃止した関数など

# 各indexを辞書として保存する
function _pack(allinds, pinds, ginds)
  dinds = sort!(setdiff(allinds, ginds))
  pinds = !isempty(pinds) ? sort!(vcat(pinds...)) : pinds
  Dict(:data=>dinds, :pilot=>pinds, :guard=>ginds, :dummy=>Int[])
end