
#========================　試験的な関数 ==========================#

"""
    get_indices_set(frameinfo, frame)

data, pilot, nullのインデックスをDictで返す
"""
function get_indices(frameinfo, frame)
    rng = 1:size(frame,2)
    frameinfo = frameinfo[:,rng,:]
    # data, pilot, nullそれぞれのCartesian-index
    @views data_indices = collect(CartesianIndices(frameinfo[:,:,1])[frameinfo[:,:,1] .< 0])
    pilot.inds = CartesianIndices(frameinfo)[frameinfo .> 0]
    null_indices = CartesianIndices(frameinfo)[frameinfo .== 0]
    return Dict("data"=>data_indices, "pilot"=>pilot.inds, "null"=>null_indices)
end


"""
    get_data_indices(frame_index::AbstractArray{<:Integer}, idx; index_type)

送信フレームの有効なデータシンボルインデックスを返す(ただし空間次元は無視する)
"""
function get_data_indices(frame_index::AbstractArray{<:Integer}, rng=1:size(frame_index,2); indexing=:Linear)
    frame_view = view(frame_index, :, rng, 1) # 配列viewの作成
    if indexing==:Linear
        data_inds = LinearIndices(frame_view)[frame_view .< 0]
    elseif indexing==:Cartesian
        data_inds = CartesianIndices(frame_view)[frame_view .< 0]
    end
    data_inds
end

# フレームの指定したシンボルのパイロットインデックス取得
function get_pilot_indices(frame_index::AbstractArray{<:Integer}, idx; index_type)
    frame_view = view(frame_index, :, idx, :) # 配列viewの作成
    if index_type==:Cartesian
        pilot_inds = CartesianIndices(frame_view)[frame_view .> 0]
    elseif index_type==:Linear
        pilot_inds = LinearIndices(frame_view)[frame_view .> 0]
    end
    pilot_inds
end




# 事前に割り当てるデータindexを決定する
function pre_allocation()
end
