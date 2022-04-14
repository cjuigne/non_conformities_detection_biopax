#! /usr/bin/env python3

import json
 


def getPrefixedIdentifier(uri, prefixesDict={}):
	for (k,v) in prefixesDict.items():
		if uri.startswith(v):
			return uri.replace(v, k+":")
	return url



def getFullURI(curie, prefixesDict={}):
	for (k,v) in prefixesDict.items():
		if curie.startswith(k+":"):
			return curie.replace(k+":",v)
	return curie



def getLocalName(uri, prefixesDict={}):
	# FIXME: no need for prefixesDict
	#for (k,v) in prefixesDict.items():
	#	if url.startswith(v):
	#		return url.replace(v, "")
	#return url
	sep = max(0, uri.rfind('/')+1, uri.rfind('#')+1, uri.rfind(':')+1)
	return uri[sep:]



def getBaseURI(uri):
	sep = max(0, uri.rfind('/')+1, uri.rfind('#')+1, uri.rfind(':')+1)
	return uri[:sep]



def readPrefixesFromFile(jsonFilePath):
	with open(jsonFilePath) as jsonFile:
		prefixesDict = json.load(jsonFile)
	return prefixesDict



def convertPrefixesDictToSPARQL(prefixesDict):
	sparqlPrefixes = ""
	for (k,v) in prefixesDict.items():
		sparqlPrefixes += "PREFIX " + k + ": <" + v + ">\n"
	return sparqlPrefixes

