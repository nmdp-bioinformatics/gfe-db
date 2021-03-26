USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as row
MERGE (:GFE {
    gfe_name: row.gfe_name, // static property
    locus: row.locus // static property
    // allele_id: row.allele_id,
    // name: row.hla_name,
    // a_name: row.a_name,
    // sequence: row.sequence,
    // length: row.length,
    // release: row.imgt_release
})
MERGE (:IMGT_HLA {
    locus: row.locus,
    name: row.hla_name
})
WITH row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (imgt:IMGT_HLA { name: row.hla_name })
MERGE (imgt)-[rel:HAS_GFE]->(gfe)
    ON CREATE SET rel.release = row.imgt_release
WITH row
MERGE (:Sequence {
    gfe_name: row.gfe_name,
    locus: row.locus,
    sequence: row.sequence,
    length: row.length
})
WITH row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (seq:Sequence { gfe_name: row.gfe_name })
MERGE (gfe)-[:HAS_SEQUENCE]->(seq);
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_features.RELEASE.csv' as row
MERGE (:Feature {
    locus: row.locus,
    gfe_name: row.gfe_name,
    rank: row.rank,
    term: row.term,
    accession: row.accession,
    sequence: row.sequence,
    length: size(row.sequence),
    hash_code: row.hash_code
})
WITH row
MATCH (gfe:GFE { gfe_name: row.gfe_name })
MATCH (f:Feature { gfe_name: row.gfe_name })
MERGE (gfe)-[:HAS_FEATURE]->(f);
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'GEN_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (:GenomicAlignment {
            gfe_name: align_row.gfe_name,
            name: align_row.hla_name,
            a_name: align_row.a_name,
            rank: align_row.rank,
            bp_sequence: align_row.bp_sequence,
            length: align_row.length
            })
)
WITH align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'NUC_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (:NucleotideAlignment {
            gfe_name: align_row.gfe_name,
            name: align_row.hla_name,
            rank: align_row.rank,
            bp_sequence: align_row.bp_sequence,
            length: align_row.length
            })
)
WITH align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'PROT_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (:ProteinAlignment {
            gfe_name: align_row.gfe_name,
            name: align_row.hla_name,
            rank: align_row.rank,
            aa_sequence: align_row.aa_sequence,
            length: align_row.length
            })
)
WITH align_row
MATCH (gfe:GFE { gfe_name: align_row.gfe_name })
MATCH (gen:GenomicAlignment { gfe_name: align_row.gfe_name })
MATCH (nuc:NucleotideAlignment { gfe_name: align_row.gfe_name })
MATCH (prot:ProteinAlignment { gfe_name: align_row.gfe_name }) 
MERGE (gfe)-[:HAS_ALIGNMENT]->(gen)
MERGE (gfe)-[:HAS_ALIGNMENT]->(nuc) 
MERGE (gfe)-[:HAS_ALIGNMENT]->(prot);
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_groups.RELEASE.csv' as groups_row
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'G' THEN [1] 
    ELSE [] END |
        MERGE (:G {
            locus: groups_row.locus,
            allele_id: groups_row.allele_id,
            name: groups_row.hla_name,
            a_name: groups_row.a_name,
            ard_id: groups_row.ard_id,
            ard_name: groups_row.ard_name
        })
)
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lg' THEN [1] 
    ELSE [] END |
        MERGE (:lg {
            locus: groups_row.locus,
            allele_id: groups_row.allele_id,
            name: groups_row.hla_name,
            a_name: groups_row.a_name,
            ard_id: groups_row.ard_id,
            ard_name: groups_row.ard_name
        })
)
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lgx' THEN [1] 
    ELSE [] END |
        MERGE (:lgx {
            locus: groups_row.locus,
            allele_id: groups_row.allele_id,
            name: groups_row.hla_name,
            a_name: groups_row.a_name,
            ard_id: groups_row.ard_id,
            ard_name: groups_row.ard_name
        })
)
WITH groups_row
MATCH (hla:IMGT_HLA { name: groups_row.hla_name })
MATCH (_g:G { name: groups_row.hla_name }) 
MATCH (_lg:lg { name: groups_row.hla_name }) 
MATCH (_lgx:lgx { name: groups_row.hla_name })
MERGE (hla)-[:G]->(_g)
MERGE (hla)-[:lg]->(_lg)
MERGE (hla)-[:lgx]->(_lgx);
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_cds.RELEASE.csv' as cds_row
MERGE (:CDS {
    gfe_name: cds_row.gfe_name,
    bp_sequence: cds_row.bp_sequence,
    bp_length: size(cds_row.bp_sequence),
    aa_sequence: cds_row.aa_sequence,
    aa_length: size(cds_row.aa_sequence)
})
WITH cds_row
MATCH (seq:Sequence { gfe_name: cds_row.gfe_name })
MATCH (cds:CDS { gfe_name: cds_row.gfe_name })
MERGE (seq)-[:HAS_CDS]->(cds);
