using Comm
include("./SphereDetector.jl")
SD = SphereDetector

Ntx, Nrx = 4, 4
SNR_dB = 20
N0 = 10^(-SNR_dB/10)

x = 2.0 * rand(0:1, Ntx) .- 1.0
H = randn(Nrx, Ntx)
y = H * x + sqrt(N0) * randn(Nrx)
refs = [-1.0, 1.0]
qam = Qam(2)
slicer = Comm.DigitalModulation.Slicer(qam)
metrics, vlist = SD.sphere_detector(y, H, refs, slicer, output=:list)

@enter SD.sphere_detector(y, H, refs, slicer)