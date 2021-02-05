MATCH (n) DETACH DELETE n;

// NODES 
LOAD CSV WITH HEADERS // Passenger
FROM "file:///titanic_clean.csv" AS row
CREATE (p:Passenger {
	name: row.name,
    age: toFloat(row.age),
    embarked: row.embarked,
    destination: row.`home.dest`,
    home_country: row.`home.country`,
    pclass: toInteger(row.pclass),
	fare: toFloat(row.fare),
    ticket: row.ticket,
    sibsp: row.sibsp,
    parch: row.parch,
    family_size: toInteger(row.`family.size`),
    surname: row.surname,
    cabin: row.cabin,
    deck: row.deck,
    sex: row.sex,
    survived: toInteger(row.survived),
    lifeboat_no: row.boat,
    body: row.body
    });
MATCH (p:Passenger) // Class
WITH DISTINCT p.pclass as class
CREATE (c:Class { pclass: class });
MATCH (p:Passenger) // Cabin
WITH DISTINCT split(p.cabin, ' ') as cabins
UNWIND cabins AS cabin
CREATE (c:Cabin { cabin: cabin });
MATCH (p:Passenger) // Embarked
WITH DISTINCT p.embarked as embarked
CREATE (:Embarked { embarked: embarked });
MATCH (p:Passenger) // Lifeboat
WITH split(p.lifeboat_no, ' ') AS lifeboats
UNWIND lifeboats AS lifeboat
WITH DISTINCT lifeboat AS boat
CREATE (:Lifeboat { lifeboat_no: boat });
MATCH (p:Passenger) // Ticket
WITH DISTINCT p.ticket AS ticket
CREATE (:Ticket { ticket: ticket });
MATCH (p:Passenger) // Deck
WITH DISTINCT p.deck AS deck
CREATE (:Deck { deck: deck });
MATCH (p:Passenger) // Destination
WITH split(p.home_country, ', ') AS countries
UNWIND countries AS country
WITH DISTINCT country AS name
CREATE (:Country { name: name });

// RELATIONSHIPS
MATCH (p:Passenger) // IN_CLASS
WITH p, p.pclass as class
MATCH (c:Class { pclass: class })
MERGE (p)-[:IN_CLASS]->(c);
MATCH (p:Passenger) // SLEPT_IN
WITH p, split(p.cabin, ' ') as cabins
UNWIND cabins AS cabin
MATCH (c:Cabin { cabin: cabin })
MERGE (p)-[:SLEPT_IN]->(c);
MATCH (p:Passenger) // EMBARKED_FROM
WITH p, p.embarked as embarked
MATCH (e:Embarked { embarked: embarked })
MERGE (p)-[:EMBARKED_FROM]->(e);
MATCH (p:Passenger) // BOARDED
WITH p, p.lifeboat_no as lifeboat
MATCH (l:Lifeboat { lifeboat_no: lifeboat })
MERGE (p)-[:BOARDED]->(l);
MATCH (p:Passenger) // ON_DECK
WITH p, p.deck as deck
MATCH (d:Deck { deck: deck })
MERGE (p)-[:ON_DECK]->(d);
MATCH (p:Passenger) // HELD_TICKET
WITH p, p.ticket as ticket
MATCH (t:Ticket { ticket: ticket })
MERGE (p)-[:HELD_TICKET]->(t);
MATCH (p1:Passenger) // RELATED_TO
WITH p1, trim(p1.surname) as surname
MATCH (p2:Passenger { surname: surname })
WHERE p1.name <> p2.name
MERGE (p1)<-[:RELATED_TO]->(p2);
MATCH (p:Passenger) // TRAVELING_TO
WITH p, split(p.home_country, ', ') as countries
UNWIND countries AS country
MATCH (c:Country { name: country })
MERGE (p)-[:TRAVELING_TO]->(c);
MATCH (l:Lifeboat) // Rename lifeboat labels
SET l.lifeboat_no = "Lifeboat " + l.lifeboat_no;