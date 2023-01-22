
function signal_detection(y, H, N0=nothing; mod_type, dims=1, method, output_type=:hard)
    # アンテナ次元
    if dims==1
        nsym, nrx, ntx = size(y,1), size(H,2), size(H,3)
    elseif dims == 2
        nsym, nrx, ntx = size(y,1), size(H,1), size(H,3)
    end

    # 推定シンボル配列，等価雑音配列の確保
    x = Array{eltype(y)}(undef, ntx, nsym)
    if output_type==:soft
        nvars = similar(x, Float64)
    end

    # SISO detection
    if ntx==1 && ntx == 1
        output_type==:hard && return siso_zf(y, H)
        output_type==:soft && return siso_zf_soft(y, H, N0)

    # SIMO detection
    elseif ntx==1 && ntx>1
        if method==:SC
        elseif method==:EGC
        elseif method==:MRC
            output_type==:hard && return simo_mrc(y, H, N0)
        end

    # MIMO detection
    elseif ntx>1 && ntx>1
        if method==:ZF
            output_type==:hard && return routine!(MIMO.zf!, x, y, H)
            output_type==:soft && return routine!(MIMO.zf_soft!, x, nvars, y, H, N0)

        elseif method==:MMSE
            if output_type == :hard
                return routine!(MIMO.mmse!, x, y, H, N0, Dim(dims))
            else
                return routine!(MIMO.mmse_soft!, x, nvars, y, H, N0, Dim(dims))
            end

        elseif method==:VBLAST
            if args.detector == "ZF"
                N0 = nothing
                params = (slicer=x->slice(mod_type, x), detector=args.detector)
            else
                params = (slicer=x->slice(mod_type, x), detector=args.detector)
            end
            output_type==:hard && return routine!(MIMO.vblast!, x, y, H, N0; params...)
            output_type==:soft && return routine!(MIMO.vblast!, x, y, H, N0; params...)

        elseif method==:SIC
            args = merge(args, (slicer = x->slicer(mod_type, x),))
            output_type==:hard && return mimo_sic(y, H, N0, mod_type)
            output_type==:soft && return mimo_sic_soft(y, H, option)

        elseif method==:PIC

        elseif method==:MLD
            args = make_refs(mod_type)
            output_type==:hard && return routine!(MIMO.mld!, x, y, H, args...)
            output_type==:soft && return routine!(MIMO.mld!, x, y, H, args...)

        elseif method==:QRMLD
            args = merge(args, (slicer = x->slicer(mod_type, x), refs=getrefs(mod_type),))
            output_type==:hard && return routine!(qrmld!, x, y, H; args...)
            output_type==:soft && return mimo_qrmld_soft(y, H, mod_type, M)

        elseif method==:SphereDecoder
            output_type==:hard && return mimo_sphere_detection(y, H, output_type; option...)
            output_type==:soft && return mimo_sphere_detection(y, H, output_type; option...)
        else
            error("該当するmethod=$(method)が見つかりませんでした．")
        end
    end
end
