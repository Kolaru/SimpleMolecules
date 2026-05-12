using SimpleMolecules

formic = read("example/formic_acid.xyz", CartesianMolecule)
iformic = InternalCoordinateMolecule(formic)
formic2 = CartesianMolecule(iformic)
iformic2 = InternalCoordinateMolecule(formic2)
formic3 = CartesianMolecule(iformic2)

dithiane = read("example/dithiane.xyz", CartesianMolecule ; bond_tolerance = 0.3)
idithiane = InternalCoordinateMolecule(dithiane ; verbose = true) ;
dithiane2 = CartesianMolecule(idithiane)