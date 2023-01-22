module LogSumProduct

include("utils.jl")
using SparseArrays

export log_sum_product!

abstract type Check end
abstract type Var end

struct CheckNode
    id::Int64
    neighbors::Vector{Int64}
    container::SparseVector{Int64, Float64}
end

struct VarNode
    id::Int64
    value::Vector{Float64}
    neighbors::Vector{Int64}
    others::SparseVector{Int64, Float64}
end

const CheckNodes = Vector{CheckNode}
const VarNodes = Vector{VarNode}

struct Graph
    checknodes::CheckNodes
    varnodes::VarNodes
    port::SparseVector{Int64, Float64}
end

function init(graph::Graph, Lpri)
    varnodes = graph.varnodes
    # Threads.@threads 
    for varnode in varnodes
        varnode.port[varnode.neighbors] .= 0.0
        varnode.value = tanh(Lpri[varnode.id]/2.0)
    end
end

function message_passing(checknodes, varnodes)
    # var-node -> check-node
    for varnode in eachindex(varnodes) # Threads.@threads 
        for id in varnode.neighbors
            message = tanh((varnode.value - varnode.port[id])/2.0)
            checknodes[id].port[varnode.id] = 
        end
    end
    # check-node -> var-node
    for checknode in checknodes
        for id in checknode.neighbors
            # message calc 
            message = message_calc()
            # message passing
            varnodes[id].port[checknode.id] = message
        end
    end

    # finale message update
    for varnode in varnodes
        varnode.value = sum(port)
    end
end


function pass(checknodes::CheckNodes, varnodes::VarNodes)
    
end



"""
    log_sum_product!(Lapp, La, rows, columns, cnlist, bnlist, iter)
log-sum-product復号
"""
function log_sum_product!(Lpos, Lpri, rows, columns, cnlist, bnlist, iterations)
    iter = 0; cks = 1;
    while iter < iterations && cks>0
        _message_passing() # メッセージ交換
        lapp_calc!(Lapp, La, obj.rows, obj.cnlist) # 事後確率計算
        decision!(decbits, Lapps) # 判定
        cks = check_sum(decbits, obj.bnlist) # チェックサム
        iter += 1
    end
    Lapp
end



# メッセージパッシング
function _message_passing()
    if iter==0
        # _bit_to_check!(n, La[n], cnlist[n], columns[n])
        Threads.@threads for n in eachindex(cnlist)
            mes = tanh(La[n]/2.0) # sending message
            for id in neighbors
                checknodes[id][] = mes
            end
        end
    else
        Threads.@threads for n in eachindex(cnlist)
            neighbors = cnlist[n]
            Lapp = La[n]
            for id in neighbors # チェックノード集合
                Lapp += rows[id][n] # 外部メッセージ
            end
            for id in neighbors
                Lext = Lapp - rows[id][n]
                columns[n][id] = tanh(Lext/2.0)
                #mes = Lsum - rows[i][n]
                #columns[n][i] = mysign(mes)*phi(mes)
            end
        end
    end
    # checknode to bitnode
    Threads.@threads for m in eachindex(bnlist) # チェックノード
        # _check_to_bit!(m, bnlist[m], rows[m], columns)
        mes = 1.0 # 初期メッセージ
        neighbors = bnlist[m]

        for id in neighbors
            mes *= columns[id][m]
            #s *= mysign(columns[j][m])
            #mes += abs(columns[j][m])
        end
        for id in neighbors
            rows[m][id] = 2.0 * myatanh(mes / columns[id][m])
            #rows[m][i] = s*mysign(columns[i][m])*phi(mes - abs(columns[i][m]))
            if isnan(rows[m][id])
                error("A Nan value is detected!")
            end
        end
    end
end

# メッセージ伝搬: ビットノード -> チェックノード
# 初期化
function _bit_to_check!(n, La, neighbors, column)
    # mes = mysign(La[n])*phi(La[n])
    mes = tanh(La/2.0) # sending message
    for id in neighbors
        # columns[n][id] = tanh(mes/2.0)
        column[id] = mes
    end
end

# 更新
function _bit_to_check!(n, La, neighbors, column, rows)
    Lapp = La
    for id in neighbors # チェックノード集合
        Lapp += rows[id][n] # 外部メッセージ
    end
    for id in neighbors
        Lext = Lapp - rows[id][n]
        column[id] = tanh(Lext/2.0)
        #mes = Lsum - rows[i][n]
        #columns[n][i] = mysign(mes)*phi(mes)
    end
end

# メッセージ伝搬: チェックノード -> ビットノード
function _check_to_bit!(m, neighbors, row, columns)
    mes = 1.0 # 初期メッセージ
    for i in neighbors
        mes *= columns[i][m]
        #s *= mysign(columns[j][m])
        #mes += abs(columns[j][m])
    end
    for i in neighbors
        row[i] = 2.0 * myatanh(mes / columns[i][m])
        #rows[m][i] = s*mysign(columns[i][m])*phi(mes - abs(columns[i][m]))
        if isnan(row[i])
            error("A Nan value is detected!")
        end
    end
end


""" sequential_decoding(dec, La, MAXITER, period=5)

直列メッセージパッシング復号
"""
function sequential_decoding(dec, La, MAXITER, period=5)
    bnlist = dec.bnlist; cnlist = dec.cnlist
    count = ones(Int64, length(bnlist)) #
    degc = [length(i) for i in bnlist] # チェックノードの次数
    a = [zeros(length(i)-1) for i in bnlist] # forward外部値
    b = [zeros(length(i)-1) for i in bnlist] # backward外部値
    Lapps = zeros(length(cnlist)) #
    decoded_bits = zeros(Bool,length(cnlist))
    iter = 1; checksum=true
    while (iter <= MAXITER) && checksum
        sequential_decoding!(La, dec.rows, a, b, degc, count, cnlist, bnlist, iter)
        if (iter-1)%period==0 || iter==MAXITER
            Lapps = _lapp_calc!(Lapps, La, dec.rows, dec.cnlist)
            decoded_bits = _decision!(decoded_bits, Lapps)
            checksum = _check_sum!(decoded_bits, dec.bnlist)
        end
        iter += 1
    end
    # @show iter
    return Lapps
end

function sequential_decoding!(La, rows, a, b, degc, count, cnlist, bnlist, iter)
    N = length(cnlist) # 変数ノード数
    if iter==1
        for n in 1:N
            mes = La[n]
            for cn in cnlist[n]
                i = count[cn]
                if i == 1
                    a[cn][i] = mes
                    count[cn] += 1
                elseif i < degc[cn]
                    a[cn][i] = boxsum(a[cn][i-1], mes)
                    count[cn] += 1
                end
            end
        end
    end
    if (iter-1)%2==0
        # backward update
        for n in reverse(1:N)
            Lsum = La[n]
            # 変数ノードnへのメッセージ計算
            for cn in cnlist[n]
                i = count[cn]
                if i == degc[cn]
                    mes = a[cn][i-1] # 外部尤度比
                elseif 1 < i
                    mes = boxsum(a[cn][i-1],b[cn][i])
                else
                    mes = b[cn][i]
                end
                rows[cn][n] = mes　# 外部尤度比更新
                Lsum += mes
            end
            # チェックノードへのメッセージ伝搬
            for cn in cnlist[n]
                i = count[cn]
                mes = Lsum - rows[cn][n] # 新規メッセージ
                if i == degc[cn]
                    b[cn][i-1] = mes
                    count[cn] -= 1
                elseif 1 < i
                    b[cn][i-1] = boxsum(b[cn][i], mes)
                    count[cn] -= 1
                end
            end
        end
    else
        # forward update
        for n in 1:N
            Lsum = La[n]
            # 変数ノードnへのメッセージ計算
            for cn in cnlist[n]
                i = count[cn]
                if i == 1
                    mes = b[cn][i]
                elseif i < degc[cn]
                    mes = boxsum(a[cn][i-1],b[cn][i])
                else
                    mes = a[cn][i-1]
                end
                rows[cn][n] = mes　# 外部尤度比更新
                Lsum += mes
            end
            # チェックノードへのメッセージ伝搬
            for cn in cnlist[n]
                i = count[cn]
                mes = Lsum - rows[cn][n] # 新規メッセージ
                if i == 1
                    a[cn][i] = mes
                    count[cn] += 1
                elseif 1 < i < degc[cn]
                    a[cn][i] = boxsum(a[cn][i-1], mes)
                    count[cn] += 1
                end
            end
        end
    end
end


end # module
