USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MERGE (gfe:GFE { gfe_name: row.gfe_name }) // static property
ON CREATE SET gfe.locus = row.locus
MERGE (seq:Sequence { gfe_name: row.gfe_name })
ON CREATE SET seq.locus = row.locus,
    seq.sequence = row.sequence,
    seq.length = row.length
MERGE (imgt:IMGT_HLA { name: row.hla_name })
ON CREATE SET imgt.locus = row.locus
MERGE (imgt)-[rel:HAS_GFE]->(gfe)
ON CREATE SET rel.release = [row.imgt_release]
ON MATCH SET rel.release = rel.release + [row.imgt_release]
MERGE (gfe)-[:HAS_SEQUENCE]->(seq);
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///all_features.RELEASE.csv' as row
MERGE (f:Feature { gfe_name: row.gfe_name, sequence: row.sequence })
ON CREATE SET f.locus = row.locus,
    f.rank = row.rank,
    f.term = row.term,
    f.accession = row.accession,
    f.length = size(row.sequence),
    f.hash_code = row.hash_code;
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///all_features.RELEASE.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (f:Feature { gfe_name: row.gfe_name })
MERGE (gfe)-[:HAS_FEATURE]->(f);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'GEN_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (gen:GenomicAlignment { gfe_name: align_row.gfe_name })
        ON CREATE SET gen.name = align_row.hla_name,
            gen.a_name = align_row.a_name,
            gen.rank = align_row.rank,
            gen.bp_sequence = align_row.bp_sequence,
            gen.length = align_row.length
);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (gen:GenomicAlignment { gfe_name: align_row.gfe_name })
MERGE (gfe)-[:HAS_ALIGNMENT]->(gen);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'NUC_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (nuc:NucleotideAlignment { gfe_name: align_row.gfe_name })
        ON CREATE SET nuc.name = align_row.hla_name,
            nuc.rank = align_row.rank,
            nuc.bp_sequence = align_row.bp_sequence,
            nuc.length = align_row.length
);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (nuc:NucleotideAlignment { gfe_name: align_row.gfe_name })
MERGE (gfe)-[:HAS_ALIGNMENT]->(nuc);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'PROT_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (prot:ProteinAlignment { gfe_name: align_row.gfe_name })
        ON CREATE SET prot.name = align_row.hla_name,
            prot.rank = align_row.rank,
            prot.aa_sequence = align_row.aa_sequence,
            prot.length = align_row.length
);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (prot:ProteinAlignment { gfe_name: align_row.gfe_name })
MERGE (gfe)-[:HAS_ALIGNMENT]->(prot);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'G' THEN [1] 
    ELSE [] END |
        MERGE (_g:G { name: groups_row.hla_name, ard_id: groups_row.ard_id})
        ON CREATE SET _g.allele_id = groups_row.allele_id,
            _g.locus = groups_row.locus,
            _g.a_name = groups_row.a_name,
            _g.ard_name = groups_row.ard_name
);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lg' THEN [1] 
    ELSE [] END |
        MERGE (_lg:lg { name: groups_row.hla_name, ard_id: groups_row.ard_id})
        ON CREATE SET _lg.allele_id = groups_row.allele_id,
            _lg.locus = groups_row.locus,
            _lg.a_name = groups_row.a_name,
            _lg.ard_name = groups_row.ard_name
);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lgx' THEN [1] 
    ELSE [] END |
        MERGE (_lgx:lgx { name: groups_row.hla_name, ard_id: groups_row.ard_id})
        ON CREATE SET _lgx.allele_id = groups_row.allele_id,
            _lgx.locus = groups_row.locus,
            _lgx.a_name = groups_row.a_name,
            _lgx.ard_name = groups_row.ard_name
);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
MATCH (hla:IMGT_HLA { name: groups_row.hla_name })
MATCH (_g:G { name: groups_row.hla_name }) 
MATCH (_lg:lg { name: groups_row.hla_name }) 
MATCH (_lgx:lgx { name: groups_row.hla_name })
MERGE (hla)-[:G]->(_g)
MERGE (hla)-[:lg]->(_lg)
MERGE (hla)-[:lgx]->(_lgx);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_cds.RELEASE.csv' as cds_row
MERGE (cds:CDS { gfe_name: cds_row.gfe_name })
ON CREATE SET cds.bp_sequence = cds_row.bp_sequence,
    cds.bp_length = size(cds_row.bp_sequence),
    cds.aa_sequence = cds_row.aa_sequence,
    cds.aa_length = size(cds_row.aa_sequence);
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_cds.RELEASE.csv' as cds_row
MATCH (seq:Sequence { gfe_name: cds_row.gfe_name })
MATCH (cds:CDS { gfe_name: cds_row.gfe_name })
MERGE (seq)-[:HAS_CDS]->(cds);
CREATE CONSTRAINT gfe_constraint IF NOT EXISTS ON (gfe:GFE) ASSERT gfe.gfe_name IS UNIQUE;
CREATE CONSTRAINT imgt_constraint IF NOT EXISTS ON (imgt:IMGT_HLA) ASSERT imgt.name IS UNIQUE;
