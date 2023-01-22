include("../DigitalModulation.jl")

# PAM Modulation test
x = smap(Ref(Pam(2)), 0:1)
@assert x == [-1.0, 1.0]

x = smap(Ref(Pam(4)), 0:3)
@assert x == [-3.0:2.0:3.0;]

x = smap(Ref(Pam(8)), 0:7)
@assert x == [-7.0:2.0:7.0;]