module LogSumProductAlter

export log_sum_product_alter!

"""
    log_sum_product_alter!(Lch, rows, columns, cnlist, bnlist, iter)
phi関数を用いたlog-sum-productアルゴリズム
"""
function log_sum_product_alter!(La, rows, columns, cnlist, bnlist, iter)
    # bitnode to checknode
    if iter>0
        for n in eachindex(cnlist) # ビットノード Threads.@threads
            neighbor = cnlist[n]
            Lsum = Lch[n]
            for j in neighbor # チェックノード集合
                Lsum += rows[j][n] # 外部メッセージ
            end
            for i in neighbor
                columns[n][i] = tanh((Lsum - rows[i][n]) / 2.0)
            end
        end
    end
    # checknode to bitnode
    for m in eachindex(bnlist) # チェックノード
        neighbor = bnlist[m]# 周辺ビットノード
        mes = 1.0
        for j in neighbor
            mes *= columns[j][m]
        end
        for i in neighbor
            rows[m][i] = 2.0 * myatanh(mes / columns[i][m])
        end
    end
end


end
