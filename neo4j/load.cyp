// Nodes
RETURN 'Creating nodes...'

// USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.3430.csv' as row
MERGE (gfe:GFE { gfe_name: row.gfe_name }) // static property
ON CREATE SET gfe.locus = row.locus;

// USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.3430.csv' as row
MERGE (seq:Sequence { gfe_name: row.gfe_name })
ON CREATE SET seq.locus = row.locus,
    seq.sequence = row.sequence,
    seq.length = row.length;

// USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.3430.csv' as row
MERGE (imgt:IMGT_HLA { name: row.hla_name })
ON CREATE SET imgt.releases = [replace(row.imgt_release, '.', '')]
ON MATCH SET imgt.releases = imgt.releases + [replace(row.imgt_release, '.', '')];

// USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///all_features.3430.csv' as row
MERGE (f:Feature { 
    locus: row.locus,
    rank: row.rank,
    term: row.term,
    accession: row.accession
    });

// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.3430.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'GEN_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (gen:GenomicAlignment { bp_sequence: align_row.bp_sequence })
        ON CREATE SET gen.label = 'GEN',
            gen.rank = align_row.rank
);

// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.3430.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'NUC_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (nuc:NucleotideAlignment { bp_sequence: align_row.bp_sequence })
        ON CREATE SET nuc.label = 'NUC', 
            nuc.rank = align_row.rank
);

// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.3430.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'PROT_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (prot:ProteinAlignment { aa_sequence: align_row.aa_sequence })
        ON CREATE SET gen.label = 'PROT', 
            prot.rank = align_row.rank
);

// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.3430.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'G' THEN [1] 
    ELSE [] END | 
        MERGE (_g:G { ard_id: groups_row.ard_id })
        ON CREATE SET _g.label = 'G'
);

// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.3430.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lg' THEN [1] 
    ELSE [] END |
        MERGE (_lg:lg { ard_id: groups_row.ard_id })
        ON CREATE SET _lg.label = 'lg'
);

// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.3430.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lgx' THEN [1] 
    ELSE [] END |
        MERGE (:lgx { ard_id: groups_row.ard_id })
        ON CREATE SET _lgx.label = 'lgx'
);

// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_cds.3430.csv' as cds_row
MERGE (cds:CDS { gfe_name: cds_row.gfe_name })
ON CREATE SET cds.bp_sequence = cds_row.bp_sequence,
    cds.bp_length = size(cds_row.bp_sequence),
    cds.aa_sequence = cds_row.aa_sequence,
    cds.aa_length = size(cds_row.aa_sequence);

// Relationships
RETURN 'Creating (:IMGT_HLA)-[:HAS_GFE]->(:GFE) ...'
// USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.3430.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (imgt:IMGT_HLA { name: row.hla_name })
MERGE (imgt)-[rel:HAS_GFE]->(gfe)
ON CREATE SET rel.releases = [replace(row.imgt_release, '.', '')]
ON MATCH SET rel.releases = rel.releases + [replace(row.imgt_release, '.', '')];

RETURN 'Creating (:GFE)-[:HAS_SEQUENCE]->(:Sequence) ...'
// USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.3430.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (seq:Sequence { sequence: row.sequence })
MERGE (gfe)-[:HAS_SEQUENCE]->(seq);

// apoc.periodic.iterate()
RETURN 'Creating (:GFE)-[:HAS_FEATURE]->(:Feature) ...'
// USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///all_features.3430.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (f:Feature { 
    locus: row.locus,
    rank: row.rank,
    term: row.term,
    accession: row.accession
    })
MERGE (gfe)-[:HAS_FEATURE]->(f);

RETURN 'Creating (:GFE)-[:HAS_ALIGNMENT]->(:GenomicAlignment) ...'
// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.3430.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (gen:GenomicAlignment { bp_sequence: align_row.bp_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(gen);

RETURN 'Creating (:GFE)-[:HAS_ALIGNMENT]->(:NucleotideAlignment) ...'
// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.3430.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (nuc:NucleotideAlignment { bp_sequence: align_row.bp_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(nuc);

RETURN 'Creating (:GFE)-[:HAS_ALIGNMENT]->(:ProteinAlignment) ...'
// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.3430.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (prot:ProteinAlignment { aa_sequence: align_row.aa_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(prot);

// apoc.periodic.iterate()
RETURN 'Creating (:IMGT_HLA)-[:G]->(:G) ...'
// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.3430.csv' as groups_row
MATCH (hla:IMGT_HLA { name: groups_row.hla_name })
MATCH (_g:G { ard_id: groups_row.ard_id }) 
MERGE (hla)-[:G]->(_g);

// apoc.periodic.iterate()
RETURN 'Creating (:IMGT_HLA)-[:lg]->(:lg) ...'
// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.3430.csv' as groups_row
MATCH (hla:IMGT_HLA { name: groups_row.hla_name })
MATCH (_lg:lg { ard_id: groups_row.ard_ide }) 
MERGE (hla)-[:lg]->(_lg);

// apoc.periodic.iterate()
RETURN 'Creating (:IMGT_HLA)-[:lgx]->(:lgx) ...'
// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.3430.csv' as groups_row
MATCH (hla:IMGT_HLA { name: groups_row.hla_name })
MATCH (_lgx:lgx { ard_id: groups_row.ard_id })
MERGE (hla)-[:lgx]->(_lgx);

// there is no way to match sequence with cds right now
RETURN 'Creating (:Sequence)-[:HAS_CDS]->(:CDS) ...'
// USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_cds.3430.csv' as cds_row
MATCH (gfe:GFE)-[:HAS_SEQUENCE]->(seq:Sequence) // Check fields in all_cds.csv for sequence
MATCH (cds:CDS { bp_sequence: cds_row.bp_sequence })
WHERE cds_row.gfe_name = gfe.gfe_name
MERGE (seq)-[:HAS_CDS]->(cds);
