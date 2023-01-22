# PEG.jl
# LDPC符合生成のためのPEG(progressive edge growth)アルゴリズムの実装
module PEG

using LinearAlgebra

# Tree型
mutable struct Tree
    root::Int64 # 根
    leaf::Dict{Symbol,Set{Int64}} # 葉
    looplen::Float64
end

# Tree型コンストラクタ
Tree(parent) = Tree(parent, Dict(:cn=>Set([]), :bn=>Set(parent)), Inf)

"""
gen_code(; nodes, degs, degc, maxlayer=Inf, algo=:determin)

    Arguments:

    Returns:
        bnlist:    ビットノードリスト
        cnlist:    チェックノードリスト
        min_cycle: グラフの最小のサイクル長
"""
function gen_code(;nodes::Tuple, degs, degc, maxlayer=Inf, algo=:determin)
    @assert length(degs) == nodes[2]
    M, N = nodes # (Check Nodes, Bit Nodes)
    allcn = [1:M;] # チェックノード集合
    degc_count = zeros(Int64, M) # degree of check nodes
    cnlist = [Int64[] for i in 1:N] # 任意のビットノードに対するチェックノード集合のリスト
    bnlist = [Int64[] for i in 1:M] # 任意のチェックノードに対するビットノード集合のリスト
    min_cycle = Inf

    for j = 1:N # ビットノード番号
        # @show j
        for k = 1:degs[j] # ビットノードの次数分だけチェックノードを選択
            layer = Set{Int64}[]
            tree = Tree(j)
            tree, layer = _unfolding(tree.root, tree, layer, cnlist, bnlist; isloop=false, l=1, maxlayer=maxlayer) # グラフを展開する
            # すべての変数ノードに"reachble"であるか？
            if length(tree.leaf[:cn]) == M
                max = length(layer) > maxlayer/2 ? Int(maxlayer/2) : length(layer) # 最大レイヤ
                candlist = Int64[]
                for i in max:length(layer)
                    union!(candlist, layer[i])
                end
                min_cycle = 2*max
            else
                candlist = setdiff(allcn, tree.leaf[:cn]) # 候補ノード集合
                min_cycle = (min_cycle > tree.looplen) ? tree.looplen : min_cycle
            end
            candlist = candlist[sortperm(degc[candlist] .- degc_count[candlist], rev=true)]
            l = 1; id = candlist[1]
            if degc[id] - degc_count[id] == 0
                id = candlist[sortperm(degc_count[candlist])][1]
            end
            #=
            for l in candlist
                if degc_count[l] < degc[l]
                    id = l
                    break
                end
                if l == candlist[end]
                    @show max, length(layer)
                    id = candlist[1]
                end
            end=#
            if degc[id] == degc_count[id]
                @show degc_count[id], degc[id]
            end
            degc_count[id] += 1
            push!(cnlist[j], id) # ノード追加
            push!(bnlist[id], j) #
        end
    end
    return bnlist, cnlist, min_cycle
end

# 部分グラフの展開
function _unfolding(parent, tree, layer, cnlist, bnlist; isloop=false, l=1, maxlayer=Inf)

    if l > maxlayer
        tree.looplen = maxlayer*2
        return
    end
    childs = Set(Int64[]) # 子ノード集合

    # 現在のレイヤの子ノードを記録する
    for pa in parent
        union!(childs, cnlist[pa])
    end
    setdiff!(childs, tree.leaf[:cn])# ただし、親ノードは除く
    if !isempty(childs)
        push!(layer, childs)
    end
    if isempty(childs) && l == 1 # 子ノードが存在する？(展開が可能であるか？)
        return tree, layer
    else
        union!(tree.leaf[:cn], childs) # 子ノードの追加
        grandchilds = Set(Int64[]) # 孫ノード(変数ノード)
        for ch in childs
            union!(grandchilds, bnlist[ch]) # 子ノード(チェックノード)の列挙
        end
        setdiff!(grandchilds, tree.leaf[:bn]) # 子ノードを除く孫ノード集合(変数ノード集合)
        if !isempty(grandchilds) # 空である
            # 現在のレイヤを保存し、次のレイヤを展開する
            union!(tree.leaf[:bn], grandchilds) # 子ノードの追加
            _unfolding(grandchilds,tree,layer,cnlist,bnlist,l=l+1,maxlayer=maxlayer) # さらにグラフを展開
        end
    end
    return tree, layer
end


# パリティ検査行列生成
function gen_parity_mtx(bnlist, cnlist)
    M = length(bnlist); N = length(cnlist)
    parity = zeros(Bool, M, N)
    for m in axes(parity,1)
        for n in Tuple(bnlist[m])
            parity[m,n] = true
        end
    end
    return parity
end

# girth4のチェック
function check_girth4(H)
    H = float(H)
    HHt = H * H'
    HHt[diagind(HHt)] .= 0
    flag = isempty(findall(x-> x .> 1,HHt))
    flag && println("Success!")
    !flag && println("Failed!")
    return
end



end
