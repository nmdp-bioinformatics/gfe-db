

def format_uri(uri):
    return "/".join(uri.replace("https://", "neo4j+s://").replace(":7473", ":7687").split("/")[:-2])