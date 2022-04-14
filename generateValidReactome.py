#! /usr/bin/env python3

import importlib
import json
import matplotlib.pyplot as plt
import os
import pandas
import pathlib
import rdflib
import rdflib.namespace
import sparqldataframe
from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import time

reactomeVersion = 79
endpointURL = "http://localhost:3030/reactome/query"
rdfFormat = "turtle"


# step 0: retrieve Reactome
#        wget https://reactome.org/download/current/biopax.zip
#        unzip biopax.zip Homo_sapiens.owl
#        echo "Reactome version: $(expr "$(grep xml:base Homo_sapiens.owl)" : '.*http:\/\/www.reactome.org\/biopax\/\([[:digit:]]*\).*')"
#
# step 1: setup SPARQL endpoint with reactome
#        # broken with fuseki-4.4.0: 
#        # temporary fix
#        ${FUSEKI_HOME}/fuseki-server --mem --update /reactome
#        # then manually load Homo_sapiens.owl
#
# step 2: export the valid complexes
#
# step 3: export the fixed invalid complexes
#
# step 4: export reactome without complexes
#        


prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX dcterms: <http://purl.org/dc/terms/>

PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>

# Homo_sapiens-20170221.owl
#PREFIX reactome: <http://www.reactome.org/biopax/59/48887#> 
#
# Homo_sapiens-20210608.owl
#PREFIX reactome: <http://www.reactome.org/biopax/77/48887#>
#
# Homo_sapiens-20220109.owl
#PREFIX reactome: <http://www.reactome.org/biopax/79/48887#>
"""




def fixInvalidComplexes():
    ##### FIX INVALID COMPLEXES
    
    query = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
PREFIX reactome: <http://www.reactome.org/biopax/"""+str(reactomeVersion)+"""/48887#>

SELECT DISTINCT ?invalidComplex

WHERE {
  ?invalidComplex rdf:type bp3:Complex .
  ?invalidComplex bp3:component ?complexPart .
  ?complexPart rdf:type bp3:Complex .
  ?complexPart bp3:component ?invalidComplexPart .
}
"""
    
    #sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes+query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    nbInvalidComplexes = len(results["results"]["bindings"])
    i = 0
    for result in results["results"]["bindings"]:
        print("{}\t{}\t{}".format(i, time.time()-startTime, result['invalidComplex']['value']))
        validGraph = rdflib.Graph()
        validGraph.bind("bp3","http://www.biopax.org/release/biopax-level3.owl#")
        complexComponents.getRepresentationBiopaxValid(endpointURL, result['invalidComplex']['value'], prefixesDict=prefixesDict, targetGraph=validGraph, rdfFormat="turtle", biopaxFilePath="")
        with open("/home/olivier/ontology/reactome/complexes/result/reactome-v" + str(reactomeVersion) + "-" + result['invalidComplex']['value'].replace("http://www.reactome.org/biopax/" + str(reactomeVersion) + "/48887#", "") + "-valid.ttl", 'w') as rdfFile:
            rdfFile.write(validGraph.serialize(format=rdfFormat).decode('UTF-8'))
        i += 1
    
    # ${JENA_HOME}/bin/riot --time --output=Turtle result/reactome-v79-Complex*.ttl > reactome-v79-complexes-invalid-fixed.ttl




def exportValidComplexes():
    ##### EXPORT VALID PART OF REACTOME
    
    query = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
PREFIX reactome: <http://www.reactome.org/biopax/"""+str(reactomeVersion)+"""/48887#>

# Retrieve the RDF graph of the original valid complexes
# file: template-extractValidComplexes-construct.rq

CONSTRUCT {
  ?validComplex ?property ?object .

  ?stoichioIdent ?propertyStoichio ?objectStoichio .

  ?xref ?propertyXref ?objectXref .
}

WHERE {
  ?validComplex rdf:type bp3:Complex .
  FILTER NOT EXISTS {
    ?validComplex bp3:component ?invalidComplexPart .
    ?invalidComplexPart rdf:type bp3:Complex .
    ?invalidComplexPart bp3:component ?invalidComplexPartPart .
  }

  ?validComplex ?property ?object .

  OPTIONAL {
    ?validComplex bp3:componentStoichiometry ?stoichioIdent .
    ?stoichioIdent ?propertyStoichio ?objectStoichio .
  }

  OPTIONAL {
    ?validComplex bp3:xref ?xref .
    ?xref ?propertyXref ?objectXref .
  }

}
"""
    
    #sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes+query)
    sparql.setReturnFormat(rdfFormat)
    results = sparql.query().convert()
    #validGraph = rdflib.Graph()
    validGraph.parse(data=results, format=rdfFormat)
    
    # EXPORT
    with open("/home/olivier/ontology/reactome/complexes/reactome-v" + str(reactomeVersion) + "-complexes-valid.ttl", 'w') as rdfFile:
        rdfFile.write(validGraph.serialize(format=rdfFormat).decode('UTF-8'))



def deleteInvalidComplexes():
    queryPath = 'queries/template-invalidComplexesComponents-delete.rq'
    sparqlQuery = pathlib.Path(queryPath).read_text().replace('$supplementaryPrefixes$',uri_utils.convertPrefixesDictToSPARQL(prefixesDict))
    print("Calling: " + endpointURL.replace("query", "update"))
    sparql = SPARQLWrapper(endpointURL.replace("query", "update"))
    sparql.setQuery(prefixes+sparqlQuery)
    sparql.method = 'POST'
    #sparql.setReturnFormat(JSON)
    sparql.query()






startTime = time.time()

pathBiopaxComplex2graph = "/home/olivier/projects/biopaxComplex2graph"
sys.path.append(pathBiopaxComplex2graph)

import uri_utils
import complexComponents

importlib.reload(uri_utils)
importlib.reload(complexComponents)

oldWD = os.getcwd()
os.chdir(pathBiopaxComplex2graph)

prefixesDict = uri_utils.readPrefixesFromFile("defaultPrefixes-release" + str(reactomeVersion) + ".json")

validGraph = rdflib.Graph()
sparql = SPARQLWrapper(endpointURL)

fixInvalidComplexes()
#exportValidComplexes()
#deleteInvalidComplexes()

os.chdir(oldWD)

endTime = time.time()
print("Duration: {}".format(endTime - startTime))




# ${JENA_HOME}/bin/riot --time --output=Turtle reactome-v79-complexes-valid.ttl reactome-v79-complexes-invalid-fixed.ttl > reactome-v79-complexes-fixed.ttl
#
# ${JENA_HOME}/bin/riot --time --output=Turtle reactome-v79-withoutComplexes.ttl reactome-v79-complexes-fixed.ttl > reactome-v79-fixed.ttl



#PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
#PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
#PREFIX owl: <http://www.w3.org/2002/07/owl#>
#PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
#PREFIX dc: <http://purl.org/dc/elements/1.1/>
#PREFIX dcterms: <http://purl.org/dc/terms/>
#PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
#
#PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
#
#PREFIX chebi: <http://purl.obolibrary.org/obo/CHEBI_>
#PREFIX reactome: <http://www.reactome.org/biopax/79/48887#>
#PREFIX uniprot: <http://purl.uniprot.org/uniprot/>
#
#
## Retrieve the direct components of a complex that are also complexes
#SELECT (COUNT(DISTINCT ?complex) AS ?nbValidComplexes)
#
#WHERE {
  # 9106 complexes
#  ?complex rdf:type bp3:Complex .
  
  # 858 black box complexes
  #FILTER NOT EXISTS {
  #  ?complex bp3:component ?complexPart .
  #}
  
  # 8248 complexes with >=1 component
  #?complex bp3:component ?complexPart .
  
  # 8248 complexes with >=1 component, none of which is a complex with components
  #?complex bp3:component ?complexPart .
  #FILTER NOT EXISTS {
  #  ?complex bp3:component ?invalidComplexPart .
  #  ?invalidComplexPart rdf:type bp3:Complex .
  #  ?invalidComplexPart bp3:component ?invalidComplexPartPart .
  #}
  
  # 9106 valid complexes
  #FILTER NOT EXISTS {
  #  ?complex bp3:component ?invalidComplexPart .
  #  ?invalidComplexPart rdf:type bp3:Complex .
  #  ?invalidComplexPart bp3:component ?invalidComplexPartPart .
  #}
  
  # 5734 complexes with >=1 component that is a complex with components
  #?complex bp3:component ?complexPart .
  #?complexPart rdf:type bp3:Complex .
  #?complexPart bp3:component ?invalidComplexPart .
#}



