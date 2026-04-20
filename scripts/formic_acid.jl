using SimpleMolecules

mol = read("example/formic_acid.xyz", CartesianMolecule)
imol = InternalCoordinateMolecule(mol)
mol2 = CartesianMolecule(imol)
imol2 = InternalCoordinateMolecule(mol2)
mol3 = CartesianMolecule(imol2)