src = pwd() * "/Mymodule/Commun/Commun.jl"
include(src)
OFDM = Commun.Mod.OFDM

ofdm = OFDM.Ofdm(nfft=64, cpsize=16, ndim=(2,2), nulls=8, pilot_config=:lte, 
                                                          pilot_interval=2, 
                                                          pilot_space=-1)
                                                    