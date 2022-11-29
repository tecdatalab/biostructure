class TMAlingResult:
  def __init__(self,
               length_of_Chain_1_p=None,
               length_of_Chain_2_p=None,
               aligned_length_p=None,
               RMSD_p=None,
               Seq_ID_p=None,
               TM_score_normalized_Chain_1=None,
               TM_score_normalized_Chain_2=None,
               ):
    self.length_of_Chain_1 = length_of_Chain_1_p
    self.length_of_Chain_2 = length_of_Chain_2_p
    self.aligned_length = aligned_length_p
    self.RMSD = RMSD_p
    self.Seq_ID = Seq_ID_p
    self.TM_score_normalized_Chain_1 = TM_score_normalized_Chain_1
    self.TM_score_normalized_Chain_2 = TM_score_normalized_Chain_2
