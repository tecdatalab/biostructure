class PDBsAlingResult:
  def __init__(self,
               pdb1Atoms,
               pdb2Atoms,
               RMSDAfterRefinement,
               NumberAlignedAtomsAfterRefinement,
               NumberRefinementCycles,
               RMSDBeforeRefinement,
               NumberAlignedAtomsBeforeRefinement,
               RawAlignmentScore,
               NumberResiduesAligned
               ):

    self.pdb1Atoms = pdb1Atoms
    self.pdb2Atoms = pdb2Atoms
    self.RMSDAfterRefinement = RMSDAfterRefinement
    self.NumberAlignedAtomsAfterRefinement = NumberAlignedAtomsAfterRefinement
    self.NumberRefinementCycles = NumberRefinementCycles
    self.RMSDBeforeRefinement = RMSDBeforeRefinement
    self.NumberAlignedAtomsBeforeRefinement = NumberAlignedAtomsBeforeRefinement
    self.RawAlignmentScore = RawAlignmentScore
    self.NumberResiduesAligned = NumberResiduesAligned

    self.PercentageAtomsAlignedBeforeRefinement = NumberAlignedAtomsBeforeRefinement / max(pdb1Atoms, pdb2Atoms)
    self.PercentageAtomsAlignedAfterRefinement = NumberAlignedAtomsAfterRefinement / max(pdb1Atoms, pdb2Atoms)
