// RETURN '(:GFE)' AS `Creating GFE nodes...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MERGE (gfe:GFE { gfe_name: row.gfe_name }) ON CREATE SET gfe.locus = row.locus
    ',
    {batchSize:5000, parallel:true});
// RETURN '(:Submitter)' AS `Creating Submitter nodes...`;
MERGE (sub:Submitter { 
    institution: 'EMBL-EBI',
    name: '<name>',
    email: '<email>'
    });
// RETURN '(:Sequence)' AS `Creating Sequence nodes...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
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
        seq.length = row.length
    ',
    {batchSize:5000, parallel:true});
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
        accession: row.accession
        })
    ',
    {batchSize:5000, parallel:true});
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MERGE (who:WHO { name: row.hla_name })
    ON CREATE SET who.gene = row.locus
    ',
    {batchSize:5000, parallel:true});
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///all_groups.RELEASE.csv" as groups_row
    RETURN groups_row
    ','
    MATCH (who:WHO { name: groups_row.hla_name })
    FOREACH(_ IN CASE 
        WHEN groups_row.ard_name = "G" THEN [1]
        ELSE [] END | SET who.G = groups_row.ard_id)
    ',
    {batchSize:5000, parallel:true});
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///all_groups.RELEASE.csv" as groups_row
    RETURN groups_row
    ','
    MATCH (who:WHO { name: groups_row.hla_name })
    FOREACH(_ IN CASE 
        WHEN groups_row.ard_name = "lg" THEN [1]
            ELSE [] END | SET who.lg = groups_row.ard_id
    )
    ',
    {batchSize:5000, parallel:true});
// RETURN '(:WHO)-[:HAS_WHO]->(:GFE)' AS `Creating relationships...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MATCH (gfe:GFE { gfe_name: row.gfe_name })
    MATCH (who:WHO { name: row.hla_name })
    MERGE (gfe)-[rel:HAS_WHO]->(who)
    ON CREATE SET rel.releases = [replace(row.imgt_release, ".", "")]
    ON MATCH SET rel.releases = rel.releases + [replace(row.imgt_release, ".", "")]
    ',
    {batchSize:5000, parallel:false});
// RETURN '(:Submitter)-[:SUBMITTED]->(:GFE)' AS `Creating relationships...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MATCH (sub:Submitter)
    MATCH (gfe:GFE { gfe_name: row.gfe_name })
    MERGE (sub)-[s:SUBMITTED]->(gfe)
    ON CREATE SET s.submit_date = date()
    ',
    {batchSize:5000, parallel:false});
// RETURN '(:GFE)-[:HAS_SEQUENCE]->(:Sequence)' AS `Creating relationships...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///gfe_sequences.RELEASE.csv" as row
    RETURN row
    ','
    MATCH (gfe:GFE { gfe_name: row.gfe_name })
    MATCH (seq:Sequence { sequence: row.sequence })
    MERGE (gfe)-[:HAS_SEQUENCE]->(seq)
    ',
    {batchSize:5000, parallel:false});
// RETURN '(:GFE)-[:HAS_FEATURE]->(:Feature)' AS `Creating relationships...`;
CALL apoc.periodic.iterate(
    '
    LOAD CSV WITH HEADERS FROM "file:///all_features.RELEASE.csv" as row
    RETURN row
    ','
    MATCH (gfe:GFE { gfe_name: row.gfe_name })
    MATCH (f:Feature { 
        locus: row.locus,
        rank: row.rank,
        term: row.term,
        accession: row.accession
        })
    MERGE (gfe)-[:HAS_FEATURE]->(f)
    ',
    {batchSize:5000, parallel:false});
