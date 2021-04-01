// This properties should be unique
CREATE CONSTRAINT gfe_constraint IF NOT EXISTS ON (gfe:GFE) ASSERT gfe.gfe_name IS UNIQUE;
CREATE CONSTRAINT imgt_hla_constraint IF NOT EXISTS ON (imgt:IMGT_HLA) ASSERT imgt.name IS UNIQUE;
// CREATE CONSTRAINT seq_constraint IF NOT EXISTS ON (seq:Sequence) ASSERT seq.seq_id IS UNIQUE; // Hashed
CREATE INDEX seq_index FOR (seq:Sequence) ON (seq.seq_id); // Constraint slows down loading

// // Alignments
// CREATE CONSTRAINT gen_constraint IF NOT EXISTS ON (gen:GenomicAlignment) ASSERT gen.seq_id IS UNIQUE; // Hashed
// CREATE CONSTRAINT nuc_constraint IF NOT EXISTS ON (nuc:NucleotideAlignment) ASSERT nuc.seq_id IS UNIQUE; // Hashed
// CREATE CONSTRAINT prot_constraint IF NOT EXISTS ON (prot:ProteinAlignment) ASSERT prot.seq_id IS UNIQUE; // Hashed
CREATE INDEX gen_index FOR (gen:GenomicAlignment) ON (gen.seq_id); // Constraint slows down loading
CREATE INDEX nuc_index FOR (nuc:NucleotideAlignment) ON (nuc.seq_id); // Constraint slows down loading
CREATE INDEX prot_index FOR (prot:ProteinAlignment) ON (prot.seq_id); // Constraint slows down loading

// Groups
// CREATE CONSTRAINT g_constraint IF NOT EXISTS ON (_g:G) ASSERT _g.ard_id IS UNIQUE;
// CREATE CONSTRAINT lg_constraint IF NOT EXISTS ON (_lg:lg) ASSERT _lg.ard_id IS UNIQUE;
// CREATE CONSTRAINT lgx_constraint IF NOT EXISTS ON (_lgx:lgx) ASSERT _lgx.ard_id IS UNIQUE;
CREATE INDEX g_index FOR (_g:G) ON (_g.ard_id);
CREATE INDEX lg_index FOR (_lg:G) ON (_lg.ard_id);
CREATE INDEX lgx_index FOR (_lgx:G) ON (_lgx.ard_id);

// CREATE CONSTRAINT cds_constraint IF NOT EXISTS ON (cds:CDS) ASSERT cds.seq_id IS UNIQUE; // Not sure if this is unique
CREATE INDEX cds_index FOR (cds:CDS) ON (cds.gfe_name, cds.seq_id);

// Has to be index, composite constraints are only available in enterprise version
CREATE INDEX feature_index FOR (f:Feature) ON (f.locus, f.rank, f.term, f.accession);
