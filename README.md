# nuclmm: model nucleotide composition and simulate sequences with Markov chains

One way to test the probability of any given feature of a genomic sequence is to determine how often that feature occurs in sequences randomly sampled from some *sequence space*.
Purely random strings of As, Cs, Gs, and Ts are a poor reflection of real DNA, whose nucleotide composition reflects the presence of genes, regulatory sequences, transposable elements, other repetitive DNA, and various other factors.
However, it is extremely difficult to simulate all of these factors, but capturing higher-order nucleotide composition in an *N*th-order Markov chain and simulating sequences randomly from this model is straightforward.

The **nuclmm** software provides command-line tools and a Python API for training *N*th-order Markov models of nucleotide composition and then generating sequences randomly from these models.

More info coming soon!
