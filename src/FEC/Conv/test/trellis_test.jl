include(pwd() * "/Mymodule/Commun/Commun.jl")
Conv = Commun.FEC.Conv

trellis = Conv.poly2trellis(poly=["171" "133"])
