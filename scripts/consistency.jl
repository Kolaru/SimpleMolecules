using SimpleMolecules

formic = read("example/formic_acid.xyz", CartesianMolecule)
iformic = InternalCoordinateMolecule(formic ; verbose = true)
formic2 = CartesianMolecule(iformic)
iformic2 = InternalCoordinateMolecule(formic2)
formic3 = CartesianMolecule(iformic2)

dithiane = read("example/dithiane.xyz", CartesianMolecule)
idithiane = InternalCoordinateMolecule(dithiane, "S1" ; verbose = true)
dithiane2 = CartesianMolecule(idithiane ; verbose = true) ;
idithiane2 = InternalCoordinateMolecule(dithiane2, "S1" ; verbose = true)
dithiane3 = CartesianMolecule(idithiane2 ; verbose = true) ;