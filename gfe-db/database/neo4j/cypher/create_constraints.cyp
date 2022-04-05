CREATE CONSTRAINT gfe_constraint IF NOT EXISTS FOR (gfe:GFE) REQUIRE gfe.gfe_name IS UNIQUE;
CREATE CONSTRAINT who_constraint IF NOT EXISTS FOR (who:WHO) REQUIRE who.name IS UNIQUE;
CREATE CONSTRAINT submitter_constraint IF NOT EXISTS FOR (sub:Submitter) REQUIRE sub.email IS UNIQUE;
CREATE INDEX seq_index IF NOT EXISTS FOR (seq:Sequence) ON (seq.seq_id);
CREATE INDEX feature_index IF NOT EXISTS FOR (f:Feature) ON (f.locus, f.rank, f.term, f.accession);