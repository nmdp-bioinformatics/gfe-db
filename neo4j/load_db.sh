# # Pseudocode
# for release in releases; do
#     sed 's/.3360./.${release}./g'
#     cat load.cyp | database


# cat neo4j/load_new.cyp | \
#     sed 's/RELEASE/${RELEASE}/g' | \
#     docker exec --interactive gfe cypher-shell -u neo4j -p gfedb

for RELEASE in 