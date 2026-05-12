using SimpleMolecules

formic = read("example/formic_acid.xyz", CartesianMolecule)
iformic = InternalCoordinateMolecule(formic)
formic2 = CartesianMolecule(iformic)
iformic2 = InternalCoordinateMolecule(formic2)
formic3 = CartesianMolecule(iformic2)

