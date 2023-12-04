# Just some notes on isomorphism in RMG:

import rmgpy.species

# If i generate 2 species that are resonance structures of each other...
sp1 = rmgpy.species.Species(smiles='C=C[O]')
sp2 = rmgpy.species.Species(smiles='[CH2]C=O')

# the isomorphism check comes back False
print('Are C=C[O] and [CH2]C=O (species object) isomorphic?')
print(sp1.is_isomorphic(sp2))
print()

# but once you generate resonance structures, that comes back True
print('How about after we generate resonance structures?')
sp1.generate_resonance_structures()
print(sp1.is_isomorphic(sp2))
print()

# doesn't matter the order- if either of them has resonance structures, it comes back True
# But you can still do an isomorphism check on the original molecule and it will distinguish between the two:
print('Are C=C[O] and [CH2]C=O (molecule object) isomorphic?')
print(sp1.molecule[0].is_isomorphic(sp2.molecule[0]))
print()

