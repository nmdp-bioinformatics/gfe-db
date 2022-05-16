CREATE CONSTRAINT gfe_constraint IF NOT EXISTS FOR (gfe:GFE) REQUIRE gfe.gfe_name IS UNIQUE;
CREATE CONSTRAINT ipd_allele_constraint IF NOT EXISTS FOR (ipd:IPD_Allele) REQUIRE ipd.name IS UNIQUE;
CREATE CONSTRAINT submitter_constraint IF NOT EXISTS FOR (sub:Submitter) REQUIRE sub.email IS UNIQUE;
CREATE CONSTRAINT ipd_acc_constraint IF NOT EXISTS FOR (acc:IPD_Accession) REQUIRE acc.name IS UNIQUE;

CREATE INDEX seq_index IF NOT EXISTS FOR (seq:Sequence) ON (seq.seq_id);
CREATE INDEX feature_index IF NOT EXISTS FOR (f:Feature) ON (f.locus, f.rank, f.term, f.accession);