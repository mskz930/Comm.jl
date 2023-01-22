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