USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as gfe_row
MERGE (:GFE {
    locus: gfe_row.locus,
    allele_id: gfe_row.allele_id,
    name: gfe_row.hla_name,
    a_name: gfe_row.a_name,
    gfe_name: gfe_row.gfe_name,
    sequence: gfe_row.sequence,
    length: gfe_row.length
})
WITH gfe_row
MERGE (:IMGT_HLA {
    locus: gfe_row.locus,
    allele_id: gfe_row.allele_id,
    name: gfe_row.hla_name
})
WITH gfe_row
MERGE (:SEQUENCE {
    locus: gfe_row.locus,
    allele_id: gfe_row.allele_id,
    name: gfe_row.hla_name,
    gfe_name: gfe_row.gfe_name,
    sequence: gfe_row.sequence,
    length: gfe_row.length
})
WITH DISTINCT replace(gfe_row.imgt_release, ".", "") AS imgt_release
MATCH (hla:IMGT_HLA)
MATCH (gfe:GFE)
WHERE hla.name = gfe.name
MERGE (hla)-[rel:HAS_GFE]->(gfe)
ON CREATE SET rel.imgt_release = [imgt_release]
ON MATCH SET rel.imgt_release = coalesce(rel.imgt_release, []) + [imgt_release];
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///gfe_sequences.RELEASE.csv' as gfe_row
WITH DISTINCT replace(gfe_row.imgt_release, ".", "") AS imgt_release
MATCH (gfe:GFE)
MATCH (seq:SEQUENCE)
WHERE gfe.sequence = seq.sequence
MERGE (gfe)-[:HAS_SEQUENCE]->(seq)
MERGE (gfe)-[rel:HAS_ALIGNMENT]->(seq)
ON CREATE SET rel.imgt_release = [imgt_release]
ON MATCH SET rel.imgt_release = coalesce(rel.imgt_release, []) + [imgt_release];
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_features.RELEASE.csv' as feature_row
MERGE (:FEATURE {
    locus: feature_row.locus,
    allele_id: feature_row.allele_id,
    name: feature_row.hla_name,
    rank: feature_row.rank,
    term: feature_row.term,
    accession: feature_row.accession,
    sequence: feature_row.sequence,
    length: size(feature_row.sequence),
    hash_code: feature_row.hash_code
})
WITH DISTINCT replace(feature_row.imgt_release, ".", "") AS imgt_release
MATCH (gfe:GFE)
MATCH (f:FEATURE)
WHERE gfe.name = f.name
MERGE (gfe)-[rel:HAS_FEATURE]->(f);
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'GEN_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (gen:GEN_ALIGN {
            name: align_row.hla_name,
            a_name: align_row.a_name,
            rank: align_row.rank,
            bp_sequence: align_row.bp_sequence,
            length: align_row.length
            })
)
WITH DISTINCT replace(align_row.imgt_release, ".", "") AS imgt_release
MATCH (gfe:GFE)
MATCH (gen:GEN_ALIGN)
WHERE gfe.name = gen.name
MERGE (gfe)-[rel:HAS_ALIGNMENT]->(gen)
ON CREATE SET rel.imgt_release = [imgt_release]
ON MATCH SET rel.imgt_release = coalesce(rel.imgt_release, []) + [imgt_release];
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'NUC_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (:NUC_ALIGN {
            name: align_row.hla_name,
            a_name: align_row.a_name,
            rank: align_row.rank,
            bp_sequence: align_row.bp_sequence,
            length: align_row.length
            })
)
WITH DISTINCT replace(align_row.imgt_release, ".", "") AS imgt_release
MATCH (gfe:GFE)
MATCH (nuc:NUC_ALIGN)
WHERE gfe.name = nuc.name
MERGE (gfe)-[rel:HAS_ALIGNMENT]->(nuc)
ON CREATE SET rel.imgt_release = [imgt_release]
ON MATCH SET rel.imgt_release = coalesce(rel.imgt_release, []) + [imgt_release]; // need to remove version if...
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_alignments.RELEASE.csv' as align_row
FOREACH(_ IN CASE 
    WHEN align_row.label = 'PROT_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (:PROT_ALIGN {
            name: align_row.hla_name,
            a_name: align_row.a_name,
            rank: align_row.rank,
            aa_sequence: align_row.aa_sequence,
            length: align_row.length
            })
)
WITH DISTINCT replace(align_row.imgt_release, ".", "") AS imgt_release
MATCH (gfe:GFE)
MATCH (prot:PROT_ALIGN)
WHERE gfe.name = prot.name
MERGE (gfe)-[rel:HAS_ALIGNMENT]->(prot)
ON CREATE SET rel.imgt_release = [imgt_release]
ON MATCH SET rel.imgt_release = coalesce(rel.imgt_release, []) + [imgt_release]; // need to remove version if...
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
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = '2nd_FIELD' THEN [1] 
    ELSE [] END |
        MERGE (sec_field_group:`2nd_FIELD` {
            locus: groups_row.locus,
            allele_id: groups_row.allele_id,
            name: groups_row.hla_name,
            a_name: groups_row.a_name,
            ard_id: groups_row.ard_id,
            ard_name: groups_row.ard_name
        })
);
MATCH (hla:IMGT_HLA)
MATCH (group:G)
WHERE hla.name = group.name
MERGE (hla)<-[rel:G]-(group);
MATCH (hla:IMGT_HLA)
MATCH (group:lg)
WHERE hla.name = group.name
MERGE (hla)<-[rel:lg]-(group);
MATCH (hla:IMGT_HLA)
MATCH (group:lgx)
WHERE hla.name = group.name
MERGE (hla)<-[rel:lgx]-(group);
MATCH (hla:IMGT_HLA)
MATCH (group:`2nd_FIELD`)
WHERE hla.name = group.name
MERGE (hla)<-[rel:`2nd_FIELD`]-(group);
USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file:///all_cds.RELEASE.csv' as cds_row
MERGE (cds:CDS {
    allele_id: cds_row.allele_id,
    name: cds_row.hla_name,
    bp_sequence: cds_row.bp_sequence,
    bp_length: size(cds_row.bp_sequence),
    aa_sequence: cds_row.aa_sequence,
    aa_length: size(cds_row.aa_sequence)
})
WITH cds
MATCH (seq:SEQUENCE)
MATCH (cds:CDS)
WHERE seq.allele_id = cds.allele_id
MERGE (seq)-[rel:HAS_CDS]->(cds);
