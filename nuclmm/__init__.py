from nuclmm.markovchain import MarkovChain


def parse_fasta(data):
    """Stolen shamelessly from http://stackoverflow.com/a/7655072/459780."""
    name, seq = None, []
    for line in data:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))
