// RETURN '(:GFE)' AS `Creating GFE nodes...`;
USING PERIODIC COMMIT 20000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MERGE (gfe:GFE { gfe_name: row.gfe_name })
ON CREATE SET gfe.locus = row.locus;
// RETURN '(:Submitter)' AS `Creating Submitter nodes...`;
MERGE (sub:Submitter { 
    institution: 'EMBL-EBI',
    name: '<name>',
    email: '<email>'
    });
// RETURN '(:Sequence)' AS `Creating Sequence nodes...`;
USING PERIODIC COMMIT 20000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MERGE (seq:Sequence { gfe_name: row.gfe_name })
ON CREATE SET 
    seq.seq_id = row.seq_id,
    seq.locus = row.locus,
    seq.sequence = row.sequence,
    seq.length = row.length
ON MATCH SET 
    seq.seq_id = row.seq_id,
    seq.locus = row.locus,
    seq.sequence = row.sequence,
    seq.length = row.length;   
// RETURN '(:Feature)' AS `Creating Feature nodes...`;
USING PERIODIC COMMIT 20000
LOAD CSV WITH HEADERS FROM 'file:///all_features.RELEASE.csv' as row
MERGE (f:Feature { 
    locus: row.locus,
    rank: row.rank,
    term: row.term,
    accession: row.accession
    });
// RETURN '(:GenomicAlignment)' AS `Creating GenomicAlignment nodes...`;
USING PERIODIC COMMIT 20000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'GEN_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (gen:GenomicAlignment { bp_sequence: align_row.bp_sequence })
        ON CREATE SET gen.label = 'GEN',
            gen.seq_id = align_row.seq_id,
            gen.rank = align_row.rank
);
// RETURN '(:NucleotideAlignment)' AS `Creating NucleotideAlignment nodes...`;
USING PERIODIC COMMIT 20000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'NUC_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (nuc:NucleotideAlignment { bp_sequence: align_row.bp_sequence })
        ON CREATE SET nuc.label = 'NUC', 
            nuc.seq_id = align_row.seq_id,
            nuc.rank = align_row.rank
);
// RETURN '(:ProteinAlignment)' AS `Creating ProteinAlignment nodes...`;
USING PERIODIC COMMIT 20000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'PROT_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (prot:ProteinAlignment { aa_sequence: align_row.aa_sequence })
        ON CREATE SET prot.label = 'PROT', 
            prot.seq_id = align_row.seq_id,
            prot.rank = align_row.rank
);
// RETURN '(:WHO)' AS `Creating WHO nodes...`;
USING PERIODIC COMMIT 20000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MERGE (who:WHO { name: row.hla_name })
ON CREATE SET who.gene = row.locus;
USING PERIODIC COMMIT 20000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
MATCH (who:WHO { name: groups_row.hla_name })
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'G' THEN [1]
    ELSE [] END | SET who.G = groups_row.ard_id
);
USING PERIODIC COMMIT 20000 
LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
MATCH (who:WHO { name: groups_row.hla_name })
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lg' THEN [1]
        ELSE [] END | SET who.lg = groups_row.ard_id
);
// RETURN '(:WHO)-[:HAS_WHO]->(:GFE)' AS `Creating relationships...`;
USING PERIODIC COMMIT 20000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (who:WHO { name: row.hla_name })
MERGE (gfe)-[rel:HAS_WHO]->(who)
ON CREATE SET rel.releases = [replace(row.imgt_release, '.', '')]
ON MATCH SET rel.releases = rel.releases + [replace(row.imgt_release, '.', '')];
// RETURN '(:Submitter)-[:SUBMITTED]->(:GFE)' AS `Creating relationships...`;
USING PERIODIC COMMIT 20000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MATCH (sub:Submitter)
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MERGE (sub)-[s:SUBMITTED]->(gfe)
ON CREATE SET s.submit_date = date();
// RETURN '(:GFE)-[:HAS_SEQUENCE]->(:Sequence)' AS `Creating relationships...`;
USING PERIODIC COMMIT 20000
LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (seq:Sequence { sequence: row.sequence })
MERGE (gfe)-[:HAS_SEQUENCE]->(seq);
// apoc.periodic.iterate()
// RETURN '(:GFE)-[:HAS_FEATURE]->(:Feature)' AS `Creating relationships...`;
USING PERIODIC COMMIT 20000
LOAD CSV WITH HEADERS FROM 'file:///all_features.RELEASE.csv' as row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (f:Feature { 
    locus: row.locus,
    rank: row.rank,
    term: row.term,
    accession: row.accession
    })
MERGE (gfe)-[:HAS_FEATURE]->(f);
// RETURN '(:GFE)-[:HAS_ALIGNMENT]->(:GenomicAlignment)' AS `Creating relationships...`;
USING PERIODIC COMMIT 20000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (gen:GenomicAlignment { bp_sequence: align_row.bp_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(gen);
// RETURN '(:GFE)-[:HAS_ALIGNMENT]->(:NucleotideAlignment)' AS `Creating relationships...`;
USING PERIODIC COMMIT 20000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (nuc:NucleotideAlignment { bp_sequence: align_row.bp_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(nuc);
// RETURN '(:GFE)-[:HAS_ALIGNMENT]->(:ProteinAlignment)' AS `Creating relationships...`;
USING PERIODIC COMMIT 20000 
LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (prot:ProteinAlignment { aa_sequence: align_row.aa_sequence })
MERGE (gfe)-[:HAS_ALIGNMENT]->(prot);
