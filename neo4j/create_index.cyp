CREATE CONSTRAINT gfe_constraint IF NOT EXISTS ON (gfe:GFE) ASSERT gfe.gfe_name IS UNIQUE;
CREATE CONSTRAINT imgt_hla_constraint IF NOT EXISTS ON (imgt:IMGT_HLA) ASSERT imgt.name IS UNIQUE;
CREATE CONSTRAINT seq_constraint IF NOT EXISTS ON (seq:Sequence) ASSERT seq.sequence IS UNIQUE;
CREATE INDEX feature_index FOR (f:Feature) ON (f.sequence, f.gfe_name);
// CREATE CONSTRAINT feature_key_constraint IF NOT EXISTS ON (f:Feature) ASSERT (f.sequence, f.gfe_name) IS NODE KEY;
