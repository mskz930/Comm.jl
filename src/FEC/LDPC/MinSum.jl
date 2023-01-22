module MinSum

export minsum!

"""
    _minsum!(La, rows, columns, cnlist, bnlist, iter)

min-sum復号
"""
function minsum!(La, rows, columns, cnlist, bnlist, iter)
    if iter == 1
        for n in eachindex(cnlist) # ビットノード Threads.@threads
            neighbor = cnlist[n]
            mes = La[n]
            for i in neighbor
                # columns[n][i] = tanh(mes/2.0)
                columns[n][i] = mes
            end
        end
    else
        Threads.@threads for n in eachindex(cnlist) # ビットノード Threads.@threads
            neighbor = cnlist[n]
            Lsum = La[n]
            for j in neighbor # チェックノード集合
                Lsum += rows[j][n] # 外部メッセージ
            end
            for i in neighbor
                # columns[n][i] = tanh((Lsum - rows[i][n]) / 2.0)
                mes = Lsum - rows[i][n]
                columns[n][i] = mes
            end
        end
    end
    # checknode to bitnode
    Threads.@threads for m in eachindex(bnlist) # チェックノード
        neighbor = bnlist[m] # 周辺ビットノード
        for i in neighbor
            mes = Inf; s=1.0
            for j in neighbor
                if i!==j
                    s *= mysign(columns[j][m])
                    mes = abs(columns[j][m]) < mes ? abs(columns[j][m]) : mes
                end
            end
            # rows[m][i] = 2.0 * myatanh(mes)
            rows[m][i] = s * mes
        end
    end
end


end # module
