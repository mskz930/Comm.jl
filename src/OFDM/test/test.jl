using Comm, Comm.OFDM

n_bits = 4096
n_dims = (2, 2)
mod_type = QamMod(4)
n_data = div(n_bits, mod_type.bps * n_dims[1])
ofdm = OfdmMod(; n_fft=128, n_gi=16, n_dims=(2,2), n_data=n_data, mod_type=mod_type, pilot=Pilot{LTE}(dt=1))

tx_bits = randbit(n_bits)
tx_data = mod(mod_type, tx_bits)
@enter OFDM.mapping(ofdm, tx_data)