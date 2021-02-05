MATCH (n) DETACH DELETE n;
LOAD CSV WITH HEADERS 
FROM 'file:///gfe_sequences.3360.csv' as gfe_row
FIELDTERMINATOR ','
// GFE nodes
MERGE (gfe:GFE {
    locus: gfe_row.locus,
    alleleId: gfe_row.alleleId,
    HLA_name: gfe_row.HLA_name,
    A_name: gfe_row.A_name,
    GFE_name: gfe_row.GFE_name,
    sequence: gfe_row.sequence,
    length: gfe_row.length
})
// SEQUENCE nodes
WITH gfe_row
MERGE (sequence:SEQUENCE {
    locus: gfe_row.locus,
    alleleId: gfe_row.alleleId,
    GFE_name: gfe_row.GFE_name,
    sequence: gfe_row.sequence,
    length: gfe_row.length
})
// FEATURE nodes
WITH max(1) AS dummy
LOAD CSV WITH HEADERS 
FROM 'file:///features.3360.csv' as feature_row
FIELDTERMINATOR ','
MERGE (feature:FEATURE {
    locus: feature_row.locus,
    alleleId: feature_row.alleleId,
    HLA_name: feature_row.HLA_name,
    rank: feature_row.rank,
    term: feature_row.term,
    accession: feature_row.accession,
    sequence: feature_row.sequence,
    length: size(feature_row.sequence),
    hash_code: feature_row.hash_code
});
// Alignments
WITH max(1) AS dummy
LOAD CSV WITH HEADERS 
FROM 'file:///alignments.3360.csv' as align_row
FIELDTERMINATOR ','
// GEN_ALIGN nodes
FOREACH(_ IN CASE 
    WHEN align_row.label = 'GEN_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (gen_align:GEN_ALIGN {
            A_name: align_row.A_name,
            rank: align_row.rank,
            sequence: align_row.sequence,
            length: align_row.length
            })
)
// NUC_ALIGN nodes
FOREACH(_ IN CASE 
    WHEN align_row.label = 'NUC_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (nuc_align:NUC_ALIGN {
            A_name: align_row.A_name,
            rank: align_row.rank,
            sequence: align_row.sequence,
            length: align_row.length
            })
)
// PROT_ALIGN nodes
FOREACH(_ IN CASE 
    WHEN align_row.label = 'PROT_ALIGN' THEN [1] 
    ELSE [] END |
        MERGE (prot_align:PROT_ALIGN {
            A_name: align_row.A_name,
            rank: align_row.rank,
            sequence: align_row.sequence,
            length: align_row.length
            })
);
// Groups
WITH max(1) AS dummy
LOAD CSV WITH HEADERS 
FROM 'file:///groups.3360.csv' as groups_row
FIELDTERMINATOR ','
// G nodes
FOREACH(_ IN CASE 
    WHEN groups_row.ARD_name = 'G' THEN [1] 
    ELSE [] END |
        MERGE (g_group:G {
            locus: groups_row.locus,
            alleleId: groups_row.alleleId,
            A_name: groups_row.A_name,
            ARD_id: groups_row.ARD_id,
            ARD_name: groups_row.ARD_name
        })
)
// lg nodes
FOREACH(_ IN CASE 
    WHEN groups_row.ARD_name = 'lg' THEN [1] 
    ELSE [] END |
        MERGE (lg_group:lg {
            locus: groups_row.locus,
            alleleId: groups_row.alleleId,
            A_name: groups_row.A_name,
            ARD_id: groups_row.ARD_id,
            ARD_name: groups_row.ARD_name
        })
)
// lgx nodes
FOREACH(_ IN CASE 
    WHEN groups_row.ARD_name = 'lgx' THEN [1] 
    ELSE [] END |
        MERGE (lgx_group:lgx {
            locus: groups_row.locus,
            alleleId: groups_row.alleleId,
            A_name: groups_row.A_name,
            ARD_id: groups_row.ARD_id,
            ARD_name: groups_row.ARD_name
        })
)
// lgx nodes
FOREACH(_ IN CASE 
    WHEN groups_row.ARD_name = '2nd_FIELD' THEN [1] 
    ELSE [] END |
        MERGE (sec_field_group:`2nd_FIELD` {
            locus: groups_row.locus,
            alleleId: groups_row.alleleId,
            A_name: groups_row.A_name,
            ARD_id: groups_row.ARD_id,
            ARD_name: groups_row.ARD_name
        })
);