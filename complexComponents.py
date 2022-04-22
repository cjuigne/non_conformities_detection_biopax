#! /usr/bin/env python3

import argparse
import collections
import distutils.util
import graphviz
import pathlib
import rdflib
import rdflib.namespace
from SPARQLWrapper import SPARQLWrapper, JSON

import uri_utils


def getComponentsDirectStoichiometry(dataSource, complexIdent, prefixesDict={}):
	"""
	Compute a list of dictionnaries of the complex direct components and their stoichiometric coefficients.

	Parameters
	----------
	dataSource : str
		Either the URL of the SPARQL endpoint or the path of the file.
	complexIdent : str
		Identifier (either CURIe or full URI) of the molecular complex.
	prefixesDict : dict of {str : str}, default {}
		Prefixes names and values.

	Returns
	-------
	listComponentStoichios : list of {'complex':???, 'complexPart':???, 'stoichioValue':???, 'stoichioIdent':???}
		List of dictionnaries of the complex components' URIs and their stoichiometric coefficients.
	"""
	queryPath = 'queries/template-complexComponentsWithStoichiometricCoefficient.rq'
	sparqlQuery = pathlib.Path(queryPath).read_text().replace('$supplementaryPrefixes$',uri_utils.convertPrefixesDictToSPARQL(prefixesDict))
	if complexIdent.startswith('http'):
		sparqlQuery = sparqlQuery.replace('$complexValue$','<' + complexIdent + '>')
	else:
		sparqlQuery = sparqlQuery.replace('$complexValue$',complexIdent)

	solutions = []
	if dataSource.startswith("http"):
		sparql = SPARQLWrapper(dataSource)
		sparql.setQuery(sparqlQuery)
		sparql.setReturnFormat(JSON)
		results = sparql.query().convert()
		for result in results["results"]["bindings"]:
			currentSolution = {}
			currentSolution['complex'] = uri_utils.getFullURI(result['complex']['value'], prefixesDict)
			currentSolution['complexPart'] = uri_utils.getFullURI(result['complexPart']['value'], prefixesDict)
			if 'stoichioValue' in result.keys():
				currentSolution['stoichioValue'] = result['stoichioValue']['value']
			else:
				currentSolution['stoichioValue'] = "1"
			if 'stoichioIdent' in result.keys():
				currentSolution['stoichioIdent'] = result['stoichioIdent']['value']
			else:
				currentSolution['stoichioIdent'] = "UnknownStoichiometry"
			solutions.append(currentSolution)
	else:
		# TODO: guess format from file extension
		g = rdflib.Graph().parse(dataSource, format="text/turtle")
		for result in g.query(sparqlQuery):
			currentSolution = {}
			currentSolution['complex'] = uri_utils.getFullURI(str(result['complex']), prefixesDict)
			currentSolution['complexPart'] = uri_utils.getFullURI(str(result['complexPart']), prefixesDict)
			currentSolution['stoichioValue'] = result['stoichioValue']
			currentSolution['stoichioIdent'] = result['stoichioIdent']
			solutions.append(currentSolution)
	return solutions


def getComponentsStoichiometry(dataSource, complexIdent, prefixesDict={}):
	"""
	Compute a dictionnary of the complex (direct and indirect) components and their stoichiometric coefficients.

	Parameters
	----------
	dataSource : str
		Either the URL of the SPARQL endpoint or the path of the file.
	complexIdent : str
		Identifier (either CURIe or full URI) of the molecular complex.
	prefixesDict : dict of {str : str}, default {}
		Prefixes names and values.

	Returns
	-------
	dictComponents : dict of {str : float}
		Dictionnary of the complex components' URIs and their stoichiometric coefficients.
	"""
	solutions = getComponentsDirectStoichiometry(dataSource, complexIdent, prefixesDict)

	if len(solutions) == 0:
		return {uri_utils.getFullURI(complexIdent, prefixesDict): 1.0}
	dictComponents = collections.defaultdict(int)
	for currentSolution in solutions:
		componentIdent = currentSolution['complexPart']
		componentStoichio = 1.0
		if 'stoichioValue' in currentSolution.keys():
			componentStoichio = float(currentSolution['stoichioValue'])
		for (k,v) in getComponentsStoichiometry(dataSource, componentIdent, prefixesDict).items():
			dictComponents[k] += componentStoichio * v
	return dictComponents




def getRepresentationBiopaxOriginal(dataSource, complexIdent, prefixesDict={}, targetGraph=None, rdfFormat="turtle", biopaxFilePath="", recursiveComplexComponents=False):
	"""
	Retrieve an RDF representation of a molecular complex consisting of all the triples having the complex as subject.

	Parameters
	----------
	dataSource : str
		Either the URL of the SPARQL endpoint or the path of the file.
	complexIdent : str
		Identifier (either CURIe or full URI) of the molecular complex.
	prefixesDict : dict of {str : str}, default {}
		Prefixes names and values.
	targetGraph : rdflib.Graph
		Graph in which to write the RDF representation. A new graph is created if None, default None.
	rdfFormat : str
		Output RDF format, default "turtle".
	biopaxFilePath : str
		File path in which to write the RDF representation if different from empty string, default "".
	recursiveComplexComponents : boolean
		recursively add to the BioPAX representation the molecular complex components that are also molecular complexes, default: False

	Returns
	-------
	string with RDF serialization.
	"""
	if targetGraph == None:
		resultGraph = rdflib.Graph()
	else:
		resultGraph = targetGraph
	queryPath = 'queries/template-complexAllRelationsAndValues-construct.rq'
	sparqlQuery = pathlib.Path(queryPath).read_text().replace('$supplementaryPrefixes$',uri_utils.convertPrefixesDictToSPARQL(prefixesDict))
	if complexIdent.startswith('http'):
		sparqlQuery = sparqlQuery.replace('$complexValue$','<' + complexIdent + '>')
	else:
		sparqlQuery = sparqlQuery.replace('$complexValue$',complexIdent)

	if dataSource.startswith("http"):
		sparql = SPARQLWrapper(dataSource)
		sparql.setQuery(sparqlQuery)
		sparql.setReturnFormat(rdfFormat)
		results = sparql.query().convert()
		resultGraph += rdflib.Graph().parse(data=results, format=rdfFormat)
		resultGraph.bind("bp3","http://www.biopax.org/release/biopax-level3.owl#")
		for (k,v) in prefixesDict.items():
			resultGraph.bind(k,v)
	else:
		# TODO: guess format from file extension
		sourceGraph = rdflib.Graph().parse(dataSource, format="text/turtle")
		results = sourceGraph.query(sparqlQuery, initNs=prefixesDict)
		resultGraph += rdflib.Graph().parse(data=results.serialize(format='xml'))
		resultGraph.bind("bp3","http://www.biopax.org/release/biopax-level3.owl#")
		for (k,v) in prefixesDict.items():
			resultGraph.bind(k,v)

	stoichios = getComponentsDirectStoichiometry(dataSource, complexIdent, prefixesDict)
	bp3NS = rdflib.Namespace("http://www.biopax.org/release/biopax-level3.owl#")
	for currentStochio in stoichios:
		resultGraph.add((rdflib.URIRef(currentStochio['stoichioIdent']), rdflib.namespace.RDF['type'], bp3NS['Stoichiometry']))
		resultGraph.add((rdflib.URIRef(currentStochio['stoichioIdent']), bp3NS['physicalEntity'], rdflib.URIRef(currentStochio['complexPart'])))
		resultGraph.add((rdflib.URIRef(currentStochio['stoichioIdent']), bp3NS['stoichiometricCoefficient'], rdflib.Literal(currentStochio['stoichioValue'], datatype=rdflib.namespace.XSD.float)))

	if recursiveComplexComponents:
		queryPath = 'queries/template-complexComponentsComplexes.rq'
		sparqlQuery = pathlib.Path(queryPath).read_text().replace('$supplementaryPrefixes$',uri_utils.convertPrefixesDictToSPARQL(prefixesDict))
		if complexIdent.startswith('http'):
			sparqlQuery = sparqlQuery.replace('$complexValue$','<' + complexIdent + '>')
		else:
			sparqlQuery = sparqlQuery.replace('$complexValue$',complexIdent)
		if dataSource.startswith("http"):
			#sparql = SPARQLWrapper(dataSource)
			sparql.setQuery(sparqlQuery)
			sparql.setReturnFormat(JSON)
			results = sparql.query().convert()
			for result in results["results"]["bindings"]:
				subComplex = result['complexPart']['value']
				getRepresentationBiopaxOriginal(dataSource, subComplex, prefixesDict, targetGraph=resultGraph, rdfFormat=rdfFormat, biopaxFilePath="", recursiveComplexComponents=recursiveComplexComponents)
		else:
			for result in sourceGraph.query(sparqlQuery, initNs=prefixesDict):
				subComplex = result['complexPart']['value']
				getRepresentationBiopaxOriginal(dataSource, subComplex, prefixesDict, targetGraph=resultGraph, rdfFormat=rdfFormat, biopaxFilePath="", recursiveComplexComponents=recursiveComplexComponents)

	if biopaxFilePath != "":
		with open(biopaxFilePath, 'w') as rdfFile:
			rdfFile.write(resultGraph.serialize(format=rdfFormat).encode().decode('UTF-8'))
	return resultGraph.serialize(format=rdfFormat).encode().decode('UTF-8')



def addToGraphBiopaxOriginal(dataSource, complexIdent, prefixesDict={}, rdfFormat="turtle", targetGraph=None, recursiveComplexComponents=False):
	"""
	Add to a graph the RDF representation of a molecular complex consisting of all the triples having the complex as subject.

	Parameters
	----------
	dataSource : str
		Either the URL of the SPARQL endpoint or the path of the file.
	complexIdent : str
		Identifier (either CURIe or full URI) of the molecular complex.
	prefixesDict : dict of {str : str}, default {}
		Prefixes names and values.
	rdfFormat : str
		Output RDF format, default "turtle".
	targetGraph : rdflib.Graph
		Graph in which to write the RDF representation. A new graph is created if None, default None.

	Returns
	-------
	string with RDF serialization.
	"""
	queryPath = 'queries/template-complexAllRelationsAndValues-construct.rq'
	sparqlQuery = pathlib.Path(queryPath).read_text().replace('$supplementaryPrefixes$',uri_utils.convertPrefixesDictToSPARQL(prefixesDict))
	if complexIdent.startswith('http'):
		sparqlQuery = sparqlQuery.replace('$complexValue$','<' + complexIdent + '>')
	else:
		sparqlQuery = sparqlQuery.replace('$complexValue$',complexIdent)

	if targetGraph == None:
		resultGraph = rdflib.Graph()
	else:
		resultGraph = targetGraph
	if dataSource.startswith("http"):
		sparql = SPARQLWrapper(dataSource)
		sparql.setQuery(sparqlQuery)
		sparql.setReturnFormat(rdfFormat)
		results = sparql.query().convert()
		resultGraph.parse(data=results, format=rdfFormat)
	else:
		# TODO: guess format from file extension
		sourceGraph = rdflib.Graph().parse(dataSource, format="text/turtle")
		results = sourceGraph.query(sparqlQuery, initNs=prefixesDict)
		resultGraph.parse(data=results.serialize(format='xml'))
		resultGraph.bind("bp3","http://www.biopax.org/release/biopax-level3.owl#")
		for (k,v) in prefixesDict.items():
			resultGraph.bind(k,v)	
	return resultGraph


def getRepresentationBiopaxValid(dataSource, complexIdent, prefixesDict={}, targetGraph=None, rdfFormat="turtle", biopaxFilePath=""):
	"""
	Generate a BioPAX-compliant RDF representation of a molecular complex such that none of its components that is also a complex has components of its own (i.e. so that all its complex components are black box complexes, cf. BioPAX specification p48).

	Parameters
	----------
	dataSource : str
		Either the URL of the SPARQL endpoint or the path of the file.
	complexIdent : str
		Identifier (either CURIe or full URI) of the molecular complex.
	prefixesDict : dict of {str : str}, default {}
		Prefixes names and values.
	targetGraph : rdflib.Graph
		Graph in which to write the RDF representation. A new graph is created if None, default None.
	rdfFormat : str
		Output RDF format, default "turtle".
	biopaxFilePath : str
		File path in which to write the RDF representation if different from empty string, default "".

	Returns
	-------
	string with RDF serialization.
	"""
	if targetGraph == None:
		resultGraph = rdflib.Graph()
	else:
		resultGraph = targetGraph
	queryPath = 'queries/template-complexAllRelationsAndValuesExceptComponents-construct.rq'
	sparqlQuery = pathlib.Path(queryPath).read_text().replace('$supplementaryPrefixes$',uri_utils.convertPrefixesDictToSPARQL(prefixesDict))
	if complexIdent.startswith('http'):
		sparqlQuery = sparqlQuery.replace('$complexValue$','<' + complexIdent + '>')
	else:
		sparqlQuery = sparqlQuery.replace('$complexValue$',complexIdent)

	if dataSource.startswith("http"):
		sparql = SPARQLWrapper(dataSource)
		sparql.setQuery(sparqlQuery)
		sparql.setReturnFormat(rdfFormat)
		results = sparql.query().convert()
		resultGraph += rdflib.Graph().parse(data=results, format=rdfFormat)
	else:
		# TODO: guess format from file extension
		sourceGraph = rdflib.Graph().parse(dataSource, format="text/turtle")
		results = sourceGraph.query(sparqlQuery, initNs=prefixesDict)
		resultGraph += rdflib.Graph().parse(data=results.serialize(format='xml'))
		resultGraph.bind("bp3","http://www.biopax.org/release/biopax-level3.owl#")
		for (k,v) in prefixesDict.items():
			resultGraph.bind(k,v)

	bp3NS = rdflib.Namespace("http://www.biopax.org/release/biopax-level3.owl#")
	complexIdentURI = rdflib.URIRef(uri_utils.getFullURI(complexIdent, prefixesDict))
	if '' in prefixesDict.keys():
		defaultNS = rdflib.Namespace(prefixesDict[''])
	else:
		defaultNS = rdflib.Namespace(uri_utils.getBaseURI(complexIdent))
	for (currentComponent, currentStoichiometricCoeff) in getComponentsStoichiometry(dataSource, complexIdent, prefixesDict).items():
		if currentComponent == complexIdent:
			continue
		getRepresentationBiopaxValid(dataSource, currentComponent, prefixesDict, resultGraph, rdfFormat="turtle", biopaxFilePath="")
		resultGraph.add((complexIdentURI, bp3NS['component'], rdflib.URIRef(currentComponent)))
		stoichiometryNodeID = 'Stoichiometry-' + uri_utils.getLocalName(complexIdent, prefixesDict) + '-' + uri_utils.getLocalName(currentComponent, prefixesDict)
		resultGraph.add((complexIdentURI, bp3NS['componentStoichiometry'], defaultNS[stoichiometryNodeID]))
		resultGraph.add((defaultNS[stoichiometryNodeID], rdflib.namespace.RDF['type'], bp3NS['Stoichiometry']))
		resultGraph.add((defaultNS[stoichiometryNodeID], bp3NS['physicalEntity'], rdflib.URIRef(currentComponent)))
		resultGraph.add((defaultNS[stoichiometryNodeID], bp3NS['stoichiometricCoefficient'], rdflib.Literal(currentStoichiometricCoeff, datatype=rdflib.namespace.XSD.float)))
	if biopaxFilePath != "":
		with open(biopaxFilePath, 'w') as rdfFile:
			rdfFile.write(resultGraph.serialize(format=rdfFormat).encode().decode('UTF-8'))
	return(resultGraph.serialize(format=rdfFormat).encode().decode('UTF-8'))




def main():
	parser = argparse.ArgumentParser(description='Generate a BioPAX-compliant representation of a BioPAX-level3 molecular complex')
	parser.add_argument("dataSource", help="either URL (starting with 'http') of the SPARQL endpoint having molecular complexes in BioPAX-level3 format, or file path (starting with anything but 'http') of the RDF file")	# 'http://localhost:3030/reactome/query'
	parser.add_argument("complexIdentifier", help="identifier of the molecular complex (CURIE or full URI)")	# 'reactome:Complex1'
	parser.add_argument("jsonPrefixesFilePath", help="path for the json file containing a dictionnary which keys are prefixes and associated values are URIs (default: empty)", nargs='?', default="")	# defaultPrefixes-release77.json
	parser.add_argument("dotFilePath", help="path for the .dot file to be generated (default: guess from complexIdentifier)", nargs='?', default="")
	parser.add_argument('--recursiveComplexComponents', help='recursively add to the BioPAX representation the molecular complex components that are also molecular complexes (default: False)', nargs='?', default=False, choices=('True','False'))
	args = parser.parse_args()
	prefixesDict = {}
	if args.jsonPrefixesFilePath != "":
		prefixesDict = uri_utils.readPrefixesFromFile(args.jsonPrefixesFilePath)


	getRepresentationBiopaxOriginal(args.dataSource, args.complexIdentifier, prefixesDict=prefixesDict, targetGraph=None, rdfFormat="turtle", biopaxFilePath=args.complexIdentifier.replace("reactome:", "")+"-orig.ttl", recursiveComplexComponents=bool(distutils.util.strtobool(str(args.recursiveComplexComponents))))

	getRepresentationBiopaxValid(args.dataSource, args.complexIdentifier, prefixesDict=prefixesDict, targetGraph=None, rdfFormat="turtle", biopaxFilePath=args.complexIdentifier.replace("reactome:", "")+"-valid.ttl")



if __name__ == "__main__":
	main()

