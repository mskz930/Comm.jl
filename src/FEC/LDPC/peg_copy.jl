module PEG

# 構造体
mutable struct Tree
    root::Integer # 根
    leaf::Dict{Symbol,Set} # 葉
    looplen::Number
end

# 外部コンストラクタ
function Tree(parent)
    Tree(parent, Dict(:cn=>Set([]), :bn=>Set(parent)), Inf)
end

# PEGアルゴリズム
function PEG_algorithm(nodes::Tuple, ds, maxlayer=Inf; alg_type=:determin)
    @assert length(ds) == nodes[2]
    M, N = nodes # (Check Nodes, Bit Nodes)
    allcn = collect(1:M)
    degc = zeros(Int64, M) # degree of check nodes
    cnlist = [Int64[] for i in 1:N] # 任意のビットノードに対するチェックノード集合のリスト
    bnlist = [Int64[] for i in 1:M] # 任意のチェックノードに対するビットノード集合のリスト
    min_cycle = Inf
    for j in 1:N # ビットノード番号
        for k in 1:ds[j] # ビットノードの次数分だけチェックノードを選択
            layer = []
            tree = gen_tree(j)
            tree, layer = unfolding(tree.root, tree, layer, cnlist, bnlist; isloop=false, l=1, maxlayer=maxlayer) # グラフを展開する
            # すべての変数ノードに"reachble"であるか？
            if length(tree.leaf[:cn]) == M
                max = length(layer) # 最大レイヤ
                candlist = collect(layer[max])
                min_cycle = 2*length(layer)
            else
                candlist = setdiff(allcn, tree.leaf[:cn]) # 候補ノード集合
                min_cycle = (min_cycle > tree.looplen) ? tree.looplen : min_cycle
            end
            id = candlist[argmin(degc[candlist])]
            degc[id] += 1
            push!(cnlist[j], id) # ノード追加
            push!(bnlist[id], j) #
        end
    end
    return bnlist, cnlist, min_cycle
end

# 現在のノードをrootとしてグラフを展開する
function unfolding(parent, tree, layer, cnlist, bnlist; isloop=false, l=1, maxlayer=Inf)

    if l > maxlayer
        tree.looplen = maxlayer*2
        return
    end

    # 現在のレイヤの子ノードを記録する
    childs = Set([]) # 子ノード集合
    for pa in parent
        union!(childs, cnlist[pa])
    end
    setdiff!(childs, tree.leaf[:cn])# ただし、親ノードは除く
    if !isempty(childs)
        push!(layer, childs)
    end
    # 子ノードが存在する？(展開が可能であるか？)
    if isempty(childs) && l == 1
        return tree, layer
    else
        union!(tree.leaf[:cn], childs) # 子ノードの追加
        grandchilds = Set([]) # 孫ノード(変数ノード)
        for ch in childs
            union!(grandchilds, bnlist[ch]) # 子ノード(チェックノード)の列挙
        end
        setdiff!(grandchilds, tree.leaf[:bn]) # 子ノードを除く孫ノード集合(変数ノード集合)
        if !isempty(grandchilds) # 空である
            # 現在のレイヤを保存し、次のレイヤを展開する
            union!(tree.leaf[:bn], grandchilds) # 子ノードの追加
            unfolding(grandchilds,tree,layer,cnlist,bnlist,l=l+1,maxlayer=maxlayer) # さらにグラフを展開
        end
    end
    tree, layer
end


end # module
