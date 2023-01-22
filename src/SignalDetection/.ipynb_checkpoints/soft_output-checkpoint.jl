
# 受信アンテナのシンボルに対する尤度比を計算
function approx_llr_from_list(v, list, pm, bitset)
    for n in size(list,1) # nrx
        for k in 1:bps
            # b[k,n]: n番目受信アンテナのk番目bit
            # b[k,n]==1の中でパスメトリックが最小の点
            id1 = id0 = -1
            C = 2.0
            max_pm1 = max_pm0 = C
            for m in axes(list,2)
                x = list[n,m]
                if x in bitset
                    id1, max_pm1 = ifelse(pm[m]>max_pm1, (m,pm[m]), (id1,min_pm1))
                else
                    id0, max_pm2 = ifelse(pm[m]>max_pm1, (m,pm[m]), (id1,min_pm1))
                end
            end
            if id0 < 0
                LLRvec[n] = max_pm1 - C
            elseif id1 < 0
                LLRvec[n] = C - min_pm0
            else
                LLRvec[n] = max_pm1 - max_pm0
            end
        end
    end
end
