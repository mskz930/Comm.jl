# methods
bps(m::Int) = ndigits(m, base=2) - 1
bps(qam::QamMod{M}) where M = bps(M)