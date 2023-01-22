# coef_matrix.jlのtest
using Commun
using Commun.OFDM.CPR
using Commun.Utils: dftmtx

pdp = Uniform(interval=5, L=32)
h = randn(ComplexF64, size(pdp.coefs)) .* pdp.coefs

L = size(h,1)-1
nfft, ngi = 64, 16
F = dftmtx(nfft, normalize=false)
Hisi, Hici = CPR.gen_coef_matrix(h, nfft, ngi, (1,1), domain=:time)
CPR.isi_coef_calc(h, F, L, ngi, nfft, 0, 0)

Hisi2 = similar(Hisi)
for k = 0:nfft-1
    for m = 0:nfft-1
        Hisi2[k+1,m+1] = CPR.isi_coef_calc(h, F, L, ngi, nfft, k, m)
    end
end
