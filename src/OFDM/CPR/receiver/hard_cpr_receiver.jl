# hard_cpr_reciver.jl


# hard-CPR receive
function cpr_receiver()
    # 事前処理(チャネル推定、係数行列の推定, データのshape)
    eq_data = similar(rx_data)
    @views Y  = similar(rx_data[:,1,:])
    @views X0 = similar(Y)
    @views X1 = similar(Y)

    #
    n = 0
    while iter < n_iters

        for n in 1:framelen
            @views Y = eq_data[:,n,:]
            if (n-1 % div_size == 0) && n > 0
                # re-estimate
                channel_restimate()
            end
            cpr(Y, X1, X0, H_coefs)

            # 信号検出
            demod_sig = signal_detection(Y, H, )
            sliced_syms = slice(demod_sig, parms=mod_parms)
            if n < framelen
                _remapping!(X0, slices_syms, timeindex=j)
            end
        end

        iter += 1 # iteration回数
    end
    eq_data
end

function _pre_processing()


end
