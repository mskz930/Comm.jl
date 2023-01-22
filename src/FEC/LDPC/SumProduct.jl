




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

# Sequential Sum-Product Decoding
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
