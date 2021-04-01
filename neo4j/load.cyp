// Nodes
RETURN '(:GFE)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MERGE (gfe:GFE { gfe_name: row.gfe_name }) // static property
ON CREATE SET gfe.locus = row.locus;

RETURN '(:Sequence)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MERGE (seq:Sequence { gfe_name: row.gfe_name })
ON CREATE SET seq.seq_id = row.seq_id,
    seq.locus = row.locus,
    seq.sequence = row.sequence,
    seq.length = row.length;

RETURN '(:Feature)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///all_features.RELEASE.csv' as row
MERGE (f:Feature { 
    locus: row.locus,
    rank: row.rank,
    term: row.term,
    accession: row.accession
    });

RETURN '(:GenomicAlignment)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'GEN_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (gen:GenomicAlignment { bp_sequence: align_row.bp_sequence })
        ON CREATE SET gen.label = 'GEN',
            gen.seq_id = align_row.seq_id,
            gen.rank = align_row.rank
);

RETURN '(:NucleotideAlignment)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'NUC_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (nuc:NucleotideAlignment { bp_sequence: align_row.bp_sequence })
        ON CREATE SET nuc.label = 'NUC', 
            nuc.seq_id = align_row.seq_id,
            nuc.rank = align_row.rank
);

RETURN '(:ProteinAlignment)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'PROT_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (prot:ProteinAlignment { aa_sequence: align_row.aa_sequence })
        ON CREATE SET prot.label = 'PROT', 
            prot.seq_id = align_row.seq_id,
            prot.rank = align_row.rank
);

RETURN '(:CDS)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_cds.RELEASE.csv' as cds_row
MERGE (cds:CDS { gfe_name: cds_row.gfe_name })
ON CREATE SET cds.bp_seq_id = cds_row.bp_seq_id,
    cds.bp_sequence = cds_row.bp_sequence,
    cds.bp_length = size(cds_row.bp_sequence),
    cds.aa_seq_id = cds_row.aa_seq_id,
    cds.aa_sequence = cds_row.aa_sequence,
    cds.aa_length = size(cds_row.aa_sequence);

RETURN '(:IMGT_HLA)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MERGE (imgt:IMGT_HLA { name: row.hla_name })
ON CREATE SET imgt.releases = [replace(row.imgt_release, '.', '')]
ON MATCH SET imgt.releases = imgt.releases + [replace(row.imgt_release, '.', '')];

RETURN '(:G)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'G' THEN [1] 
    ELSE [] END | 
        MERGE (_g:G { ard_id: groups_row.ard_id })
        ON CREATE SET _g.label = 'G'
);

RETURN '(:lg)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lg' THEN [1] 
    ELSE [] END |
        MERGE (_lg:lg { ard_id: groups_row.ard_id })
        ON CREATE SET _lg.label = 'lg'
);

RETURN '(:lgx)' AS `Creating nodes...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lgx' THEN [1] 
    ELSE [] END |
        MERGE (_lgx:lgx { ard_id: groups_row.ard_id })
        ON CREATE SET _lgx.label = 'lgx'
);

// Relationships
RETURN '(:IMGT_HLA)-[:HAS_GFE]->(:GFE)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (imgt:IMGT_HLA { name: row.hla_name })
MERGE (imgt)-[rel:HAS_GFE]->(gfe)
ON CREATE SET rel.releases = [replace(row.imgt_release, '.', '')]
ON MATCH SET rel.releases = rel.releases + [replace(row.imgt_release, '.', '')];

RETURN '(:GFE)-[:HAS_SEQUENCE]->(:Sequence)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (seq:Sequence { sequence: row.sequence })
MERGE (gfe)-[:HAS_SEQUENCE]->(seq);

// apoc.periodic.iterate()
RETURN '(:GFE)-[:HAS_FEATURE]->(:Feature)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000
LOAD CSV WITH HEADERS FROM 'file:///all_features.RELEASE.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (f:Feature { 
    locus: row.locus,
    rank: row.rank,
    term: row.term,
    accession: row.accession
    })
MERGE (gfe)-[:HAS_FEATURE]->(f);

RETURN '(:GFE)-[:HAS_ALIGNMENT]->(:GenomicAlignment)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (gen:GenomicAlignment { bp_sequence: align_row.bp_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(gen);

RETURN '(:GFE)-[:HAS_ALIGNMENT]->(:NucleotideAlignment)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (nuc:NucleotideAlignment { bp_sequence: align_row.bp_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(nuc);

RETURN '(:GFE)-[:HAS_ALIGNMENT]->(:ProteinAlignment)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (prot:ProteinAlignment { aa_sequence: align_row.aa_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(prot);

// apoc.periodic.iterate()
RETURN '(:IMGT_HLA)-[:G]->(:G)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
MATCH (hla:IMGT_HLA { name: groups_row.hla_name })
MATCH (_g:G { ard_id: groups_row.ard_id }) 
MERGE (hla)-[:G]->(_g);

// apoc.periodic.iterate()
RETURN '(:IMGT_HLA)-[:lg]->(:lg)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
MATCH (hla:IMGT_HLA { name: groups_row.hla_name })
MATCH (_lg:lg { ard_id: groups_row.ard_id }) 
MERGE (hla)-[:lg]->(_lg);

// apoc.periodic.iterate()
RETURN '(:IMGT_HLA)-[:lgx]->(:lgx)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
MATCH (hla:IMGT_HLA { name: groups_row.hla_name })
MATCH (_lgx:lgx { ard_id: groups_row.ard_id })
MERGE (hla)-[:lgx]->(_lgx);

RETURN '(:Sequence)-[:HAS_CDS]->(:CDS)' AS `Creating relationships...`;
USING PERIODIC COMMIT 50000 
LOAD CSV WITH HEADERS FROM 'file:///all_cds.RELEASE.csv' as cds_row
MATCH (gfe:GFE)-[:HAS_SEQUENCE]->(seq:Sequence)
MATCH (cds:CDS { bp_sequence: cds_row.bp_sequence })
WHERE cds_row.gfe_name = gfe.gfe_name
MERGE (seq)-[:HAS_CDS]->(cds);
