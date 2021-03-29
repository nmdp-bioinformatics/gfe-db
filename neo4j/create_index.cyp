CREATE CONSTRAINT gfe_constraint IF NOT EXISTS ON (gfe:GFE) ASSERT gfe.gfe_name IS UNIQUE;
CREATE CONSTRAINT imgt_hla_constraint IF NOT EXISTS ON (imgt:IMGT_HLA) ASSERT imgt.name IS UNIQUE;
// CREATE CONSTRAINT seq_constraint IF NOT EXISTS ON (seq:Sequence) ASSERT seq.sequence IS UNIQUE; // Too large
CREATE INDEX cds_index FOR (cds:CDS) ON (cds.gfe_name);
CREATE INDEX seq_index FOR (seq:Sequence) ON (seq.gfe_name);
CREATE INDEX feature_index FOR (f:Feature) ON (f.gfe_name);
CREATE INDEX gen_align_index FOR (gen:GenomicAlignment) ON (gen.gfe_name);
CREATE INDEX nuc_align_index FOR (nuc:NucleotideAlignment) ON (nuc.gfe_name);
CREATE INDEX prot_align_index FOR (prot:ProteinAlignment) ON (prot.gfe_name);
CREATE INDEX g_index FOR (_g:G) ON (_g.hla_name, _g.ard_id);
CREATE INDEX lg_index FOR (_lg:lg) ON (_lg.hla_name, _lg.ard_id);
CREATE INDEX lgx_index FOR (_lgx:lgx) ON (_lgx.hla_name, _lgx.ard_id);
// Requires enterprise license
// CREATE CONSTRAINT feature_key_constraint IF NOT EXISTS ON (f:Feature) ASSERT (f.sequence, f.gfe_name) IS NODE KEY;
