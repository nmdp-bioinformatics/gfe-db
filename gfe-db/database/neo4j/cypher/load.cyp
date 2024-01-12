// RETURN '(:GFE)' AS `Creating GFE nodes...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MERGE (gfe:GFE { name: row.gfe_name }) ON CREATE SET gfe.locus = row.locus
    ',
    {batchSize:500, parallel:true});
// RETURN '(:Submitter)' AS `Creating Submitter nodes...`;
MERGE (sub:Submitter { 
    institution: 'IPD',
    name: 'IPD-IMGT',
    url: 'https://www.ebi.ac.uk/ipd/imgt/hla/',
    email: '<email>'
    });
// RETURN '(:Sequence)' AS `Creating Sequence nodes...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MERGE (seq:Sequence { name: row.gfe_name }) 
    ON CREATE SET 
        seq.locus = row.locus,
        seq.sequence = row.sequence,
        seq.length = row.length
    ON MATCH SET 
        seq.locus = row.locus,
        seq.sequence = row.sequence,
        seq.length = row.length
    ',
    {batchSize:500, parallel:true});
// RETURN '(:Feature)' AS `Creating Feature nodes...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///all_features.RELEASE.csv" as row
    RETURN row
    ','
    MERGE (f:Feature { 
        locus: row.locus,
        rank: row.rank,
        term: row.term,
        accession: row.accession,
        sequence: row.sequence
        })
    ',
    {batchSize:500, parallel:true});
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MERGE (ipd:IPD_Allele { name: row.hla_name })
    ON CREATE SET ipd.gene = row.locus
    MERGE (a:IPD_Accession {
        name: row.acc_name
    })
    ',
    {batchSize:500, parallel:true});
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///all_groups.RELEASE.csv" as groups_row
    RETURN groups_row
    ','
    MATCH (ipd:IPD_Allele { name: groups_row.hla_name })
    FOREACH(_ IN CASE 
        WHEN groups_row.ard_name = "G" THEN [1]
        ELSE [] END | SET ipd.G = groups_row.ard_id)
    ',
    {batchSize:500, parallel:true});
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///all_groups.RELEASE.csv" as groups_row
    RETURN groups_row
    ','
    MATCH (ipd:IPD_Allele { name: groups_row.hla_name })
    FOREACH(_ IN CASE 
        WHEN groups_row.ard_name = "lg" THEN [1]
            ELSE [] END | SET ipd.lg = groups_row.ard_id
    )
    ',
    {batchSize:500, parallel:true});
// RETURN '(:IPD_Allele)-[:HAS_IPD_ALLELE]->(:GFE)' AS `Creating relationships...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MATCH (gfe:GFE { name: row.gfe_name })
    MATCH (ipd:IPD_Allele { name: row.hla_name })
    MATCH (acc:IPD_Accession { name: row.acc_name})
    MERGE (gfe)-[rel:HAS_IPD_ALLELE]->(ipd)
    ON CREATE SET rel.releases = [replace(row.imgt_release, ".", "")]
    ON MATCH SET rel.releases = apoc.coll.toSet(rel.releases + [replace(row.imgt_release, ".", "")])
    MERGE (gfe)-[acc_rel:HAS_IPD_ACCESSION]->(acc)
    ON CREATE SET acc_rel.release = row.imgt_release
    ',
    {batchSize:500, parallel:true});
// RETURN '(:Submitter)-[:SUBMITTED]->(:GFE)' AS `Creating relationships...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MATCH (sub:Submitter)
    MATCH (gfe:GFE { name: row.gfe_name })
    MERGE (sub)-[s:SUBMITTED]->(gfe)
    ON CREATE SET s.submit_date = date()
    ',
    {batchSize:500, parallel:false});
// RETURN '(:GFE)-[:HAS_SEQUENCE]->(:Sequence)' AS `Creating relationships...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MATCH (gfe:GFE { name: row.gfe_name })
    MATCH (seq:Sequence { sequence: row.sequence })
    MERGE (gfe)-[:HAS_SEQUENCE]->(seq)
    ',
    {batchSize:500, parallel:false});
// RETURN '(:GFE)-[:HAS_FEATURE]->(:Feature)' AS `Creating relationships...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///all_features.RELEASE.csv" as row
    RETURN row
    ','
    MATCH (gfe:GFE { name: row.gfe_name })
    MATCH (f:Feature { 
        locus: row.locus,
        rank: row.rank,
        term: row.term,
        accession: row.accession
        })
    MERGE (gfe)-[:HAS_FEATURE]->(f)
    ',
    {batchSize:500, parallel:false});
