using Commun
using BenchmarkTools

G, H = LDPC.load(rows=1024, cols=2048)
encoder = LDPC.Encoder(G)
b = randbit(1024)
encoded = encoder(b)

decoder = LDPC.Decoder(H)
