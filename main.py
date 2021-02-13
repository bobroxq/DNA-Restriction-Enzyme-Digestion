from utilities import restriction_enzymes

dna_mw = {'a':313.2, 'c':289.2, 't':304.2, 'g':329.2}

def oligo_mw(sequence):
    mw = 0.0
    for base in sequence:
        mw  += dna_mw[base]
    return mw


def restriction_digest(sequence, enzyme):
    motif = restriction_enzymes[enzyme][0]
    cut_position = restriction_enzymes[enzyme][1]
    fragments = []
    found = 0
    last_cut = found
    search_from = last_cut
    while found != -1:
        found = sequence.find(motif, search_from)
        if found != -1:
            fragment = sequence[last_cut:found + cut_position]
            mwt = oligo_mw(fragment)
            fragments.append((fragment, mwt))
        else:
            fragment = sequence[last_cut:]
            mwt = oligo_mw(fragment)
            fragments.append((fragment, mwt))
        last_cut = found + cut_position
        search_from = last_cut + 1
    return fragments

sequence = 'atgcggatccccagtacgtacgggatccatacgt'
print('Molecular weight: ', oligo_mw(sequence))
print('Fragments: ', restriction_digest(sequence, 'bamH1'))