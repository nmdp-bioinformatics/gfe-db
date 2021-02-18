// TO DO: Replace dbversion parameter in CSV file path
MATCH (n) DETACH DELETE n;
USING PERIODIC COMMIT 10000
LOAD CSV WITH HEADERS 
FROM 'file:///gfe_sequences.3360.csv' as gfe_row
FIELDTERMINATOR ','
// GFE nodes
MERGE (gfe:GFE {
    locus: gfe_row.locus,
    allele_id: gfe_row.allele_id,
    hla_name: gfe_row.hla_name,
    a_name: gfe_row.a_name,
    gfe_name: gfe_row.gfe_name,
    imgt_release: gfe_row.imgt_release,
    sequence: gfe_row.sequence,
    length: gfe_row.length
})
// IMGT_HLA nodes
WITH gfe_row
MERGE (imgt:IMGT_HLA {
    locus: gfe_row.locus,
    allele_id: gfe_row.allele_id,
    hla_name: gfe_row.hla_name
})
// SEQUENCE nodes
WITH gfe_row
MERGE (sequence:SEQUENCE {
    locus: gfe_row.locus,
    allele_id: gfe_row.allele_id,
    hla_name: gfe_row.hla_name,
    gfe_name: gfe_row.gfe_name,
    imgt_release: gfe_row.imgt_release,
    sequence: gfe_row.sequence,
    length: gfe_row.length
});
// FEATURE nodes
WITH max(1) AS dummy
USING PERIODIC COMMIT 10000
LOAD CSV WITH HEADERS 
FROM 'file:///all_features.3360.csv' as feature_row
FIELDTERMINATOR ','
MERGE (feature:FEATURE {
    imgt_release: feature_row.imgt_release,
    locus: feature_row.locus,
    allele_id: feature_row.allele_id,
    hla_name: feature_row.hla_name,
    rank: feature_row.rank,
    term: feature_row.term,
    accession: feature_row.accession,
    sequence: feature_row.sequence,
    length: size(feature_row.sequence),
    hash_code: feature_row.hash_code
});
// Alignments
// GEN_ALIGN nodes
WITH max(1) AS dummy
LOAD CSV WITH HEADERS 
FROM 'file:///all_alignments.3360.csv' as align_row
FIELDTERMINATOR ','
FOREACH(_ IN CASE 
    WHEN align_row.label = 'GEN_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (gen_align:GEN_ALIGN {
            hla_name: align_row.hla_name,
            a_name: align_row.a_name,
            rank: align_row.rank,
            bp_sequence: align_row.bp_sequence,
            imgt_release: align_row.imgt_release,
            length: align_row.length
            })
)
// (:GFE)-[:HAS_ALIGNMENT]->(GEN_ALIGN)
// TO DO: set rel.accession value
WITH align_row
MATCH (gfe:GFE)
MATCH (alignment:GEN_ALIGN)
WHERE gfe.a_name = alignment.a_name AND align_row.label = "GEN_ALIGN"
MERGE (gfe)-[rel:HAS_ALIGNMENT]->(alignment)
SET rel.imgt_release = alignment.imgt_release,
    rel.accession = "0"
// (:IMGT_HLA)-[:HAS_ALIGNMENT]->(GEN_ALIGN)
WITH align_row, alignment
MATCH (imgt:IMGT_HLA)
WHERE imgt.hla_name = alignment.hla_name AND align_row.label = "GEN_ALIGN"
MERGE (imgt)-[rel:HAS_ALIGNMENT]->(alignment)
SET rel.imgt_release = alignment.imgt_release,
    rel.accession = "0";
// NUC_ALIGN nodes
WITH max(1) AS dummy
LOAD CSV WITH HEADERS 
FROM 'file:///all_alignments.3360.csv' as align_row
FIELDTERMINATOR ','
FOREACH(_ IN CASE 
    WHEN align_row.label = 'NUC_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (nuc_align:NUC_ALIGN {
            hla_name: align_row.hla_name,
            a_name: align_row.a_name,
            rank: align_row.rank,
            bp_sequence: align_row.bp_sequence,
            imgt_release: align_row.imgt_release,
            length: align_row.length
            })
)
// (:GFE)-[:HAS_ALIGNMENT]->(NUC_ALIGN)
// TO DO: set rel.accession value
WITH align_row
MATCH (gfe:GFE)
MATCH (alignment:NUC_ALIGN)
WHERE gfe.a_name = alignment.a_name AND align_row.label = "NUC_ALIGN"
MERGE (gfe)-[rel:HAS_ALIGNMENT]->(alignment)
SET rel.imgt_release = alignment.imgt_release,
    rel.accession = "0"
// (:IMGT_HLA)-[:HAS_ALIGNMENT]->(NUC_ALIGN)
WITH align_row, alignment
MATCH (imgt:IMGT_HLA)
WHERE imgt.hla_name = alignment.hla_name AND align_row.label = "NUC_ALIGN"
MERGE (imgt)-[rel:HAS_ALIGNMENT]->(alignment)
SET rel.imgt_release = alignment.imgt_release,
    rel.accession = "0";
// PROT_ALIGN nodes
WITH max(1) AS dummy
LOAD CSV WITH HEADERS 
FROM 'file:///all_alignments.3360.csv' as align_row
FIELDTERMINATOR ','
FOREACH(_ IN CASE 
    WHEN align_row.label = 'PROT_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (prot_align:PROT_ALIGN {
            hla_name: align_row.hla_name,
            a_name: align_row.a_name,
            rank: align_row.rank,
            aa_sequence: align_row.aa_sequence,
            imgt_release: align_row.imgt_release,
            length: align_row.length
            })
)
// (:GFE)-[:HAS_ALIGNMENT]->(PROT_ALIGN)
// TO DO: set rel.accession value
WITH align_row
MATCH (gfe:GFE)
MATCH (alignment:PROT_ALIGN)
WHERE gfe.a_name = alignment.a_name AND align_row.label = "PROT_ALIGN"
MERGE (gfe)-[rel:HAS_ALIGNMENT]->(alignment)
SET rel.imgt_release = alignment.imgt_release,
    rel.accession = "0"
// (:IMGT_HLA)-[:HAS_ALIGNMENT]->(PROT_ALIGN)
WITH align_row, alignment
MATCH (imgt:IMGT_HLA)
WHERE imgt.hla_name = alignment.hla_name AND align_row.label = "PROT_ALIGN"
MERGE (imgt)-[rel:HAS_ALIGNMENT]->(alignment)
SET rel.imgt_release = alignment.imgt_release,
    rel.accession = "0";
// Groups
WITH max(1) AS dummy
LOAD CSV WITH HEADERS 
FROM 'file:///all_groups.3360.csv' as groups_row
FIELDTERMINATOR ','
// G nodes
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'G' THEN [1] 
    ELSE [] END |
        MERGE (g_group:G {
            locus: groups_row.locus,
            allele_id: groups_row.allele_id,
            hla_name: groups_row.hla_name,
            a_name: groups_row.a_name,
            ard_id: groups_row.ard_id,
            ard_name: groups_row.ard_name,
            imgt_release: groups_row.imgt_release
        })
)
// lg nodes
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lg' THEN [1] 
    ELSE [] END |
        MERGE (lg_group:lg {
            locus: groups_row.locus,
            allele_id: groups_row.allele_id,
            hla_name: groups_row.hla_name,
            a_name: groups_row.a_name,
            ard_id: groups_row.ard_id,
            ard_name: groups_row.ard_name,
            imgt_release: groups_row.imgt_release
        })
)
// lgx nodes
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = 'lgx' THEN [1] 
    ELSE [] END |
        MERGE (lgx_group:lgx {
            locus: groups_row.locus,
            allele_id: groups_row.allele_id,
            hla_name: groups_row.hla_name,
            a_name: groups_row.a_name,
            ard_id: groups_row.ard_id,
            ard_name: groups_row.ard_name,
            imgt_release: groups_row.imgt_release
        })
)
// lgx nodes
FOREACH(_ IN CASE 
    WHEN groups_row.ard_name = '2nd_FIELD' THEN [1] 
    ELSE [] END |
        MERGE (sec_field_group:`2nd_FIELD` {
            locus: groups_row.locus,
            allele_id: groups_row.allele_id,
            hla_name: groups_row.hla_name,
            a_name: groups_row.a_name,
            ard_id: groups_row.ard_id,
            ard_name: groups_row.ard_name,
            imgt_release: groups_row.imgt_release
        })
);
// CDS nodes
WITH max(1) AS dummy
USING PERIODIC COMMIT 10000
LOAD CSV WITH HEADERS 
FROM 'file:///all_cds.3360.csv' as cds_row
FIELDTERMINATOR ','
MERGE (cds:CDS {
    allele_id: cds_row.allele_id,
    hla_name: cds_row.hla_name,
    imgt_release: cds_row.imgt_release,
    bp_sequence: cds_row.bp_sequence,
    bp_length: size(cds_row.bp_sequence),
    aa_sequence: cds_row.aa_sequence,
    aa_length: size(cds_row.aa_sequence)
})
// (:SEQUENCE)-[:HAS_CDS]->(:CDS)
// TO DO: validate this relationship, confirm which value to match on
WITH cds
MATCH (seq:SEQUENCE)
MATCH (cds:CDS)
WHERE seq.allele_id = cds.allele_id
MERGE (seq)-[rel:HAS_CDS]->(cds);
// (:GFE)-[:HAS_SEQUENCE]->(SEQUENCE)
MATCH (gfe:GFE)
MATCH (seq:SEQUENCE)
WHERE gfe.sequence = seq.sequence
MERGE (gfe)-[rel:HAS_SEQUENCE]->(seq)
SET rel.imgt_release = seq.imgt_release,
    rel.accession = "0";
// (:GFE)-[:HAS_ALIGNMENT]->(SEQUENCE)
MATCH (gfe:GFE)
MATCH (seq:SEQUENCE)
WHERE gfe.sequence = seq.sequence
MERGE (gfe)-[rel:HAS_ALIGNMENT]->(seq)
SET rel.imgt_release = seq.imgt_release,
    rel.accession = "0";
// (:GFE)-[:HAS_FEATURE]->(FEATURE)
MATCH (gfe:GFE)
MATCH (f:FEATURE)
WHERE gfe.hla_name = f.hla_name
MERGE (gfe)-[rel:HAS_FEATURE]->(f)
SET rel.imgt_release = f.imgt_release,
    rel.accession = f.accession;
// (:IMGT_HLA)-[:HAS_GFE]->(:GFE)
MATCH (hla:IMGT_HLA)
MATCH (gfe:GFE)
WHERE hla.hla_name = gfe.hla_name
MERGE (hla)-[rel:HAS_GFE]->(gfe)
SET rel.imgt_release = gfe.imgt_release;
// (:IMGT_HLA)-[:HAS_SEQUENCE]->(SEQUENCE)
MATCH (hla:IMGT_HLA)
MATCH (seq:SEQUENCE)
WHERE hla.allele_id = seq.allele_id
MERGE (hla)-[rel:HAS_SEQUENCE]->(seq)
SET rel.imgt_release = seq.imgt_release,
    rel.accession = "0";
// (:IMGT_HLA)-[:HAS_ALIGNMENT]->(SEQUENCE)
MATCH (hla:IMGT_HLA)
MATCH (seq:SEQUENCE)
WHERE hla.allele_id = seq.allele_id
MERGE (hla)-[rel:HAS_ALIGNMENT]->(seq)
SET rel.imgt_release = seq.imgt_release,
    rel.accession = "0";
// (:IMGT_HLA)-[:HAS_FEATURE]->(FEATURE)
MATCH (hla:IMGT_HLA)
MATCH (f:FEATURE)
WHERE hla.hla_name = f.hla_name
MERGE (hla)-[rel:HAS_FEATURE]->(f)
SET rel.imgt_release = f.imgt_release,
    rel.accession = f.accession;
// (:IMGT_HLA)-[:G]->(GROUP)
MATCH (hla:IMGT_HLA)
MATCH (group:G)
WHERE hla.hla_name = group.hla_name
MERGE (hla)<-[rel:G]-(group)
SET rel.imgt_release = group.imgt_release;
// (:IMGT_HLA)-[:lg]->(GROUP)
MATCH (hla:IMGT_HLA)
MATCH (group:lg)
WHERE hla.hla_name = group.hla_name
MERGE (hla)<-[rel:lg]-(group)
SET rel.imgt_release = group.imgt_release;
// (:IMGT_HLA)-[:lgx]->(GROUP)
MATCH (hla:IMGT_HLA)
MATCH (group:lgx)
WHERE hla.hla_name = group.hla_name
MERGE (hla)<-[rel:lgx]-(group)
SET rel.imgt_release = group.imgt_release;
// (:IMGT_HLA)-[:`2nd_FIELD`]->(GROUP)
MATCH (hla:IMGT_HLA)
MATCH (group:`2nd_FIELD`)
WHERE hla.hla_name = group.hla_name
MERGE (hla)<-[rel:`2nd_FIELD`]-(group)
SET rel.imgt_release = group.imgt_release;