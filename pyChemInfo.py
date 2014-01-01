#!/usr/bin/env python
# pyChemInfo.py -- Informations about chemical compound.
# -*- coding: utf-8 -*-

"""
pyChemInfo

A simple small script for display informations about a chemical compound. It
takes its informations from different online databases.

All gettet informations can be save local for fast information getting on next
search.

Usage (display informations about a structure):
./pyChemInfo <Name>
"""

import zlib
import sys
import time
import re
import getopt
import os
import os.path

import abc
from abc import abstractmethod

import urllib.parse
import urllib.request

from xml.etree import ElementTree as ET


AUTHOR = "Christian Krippendorf"
EMAIL = "Kontakt@Christian-Krippendorf.de"
VERSION = "1.0"
LICENCE = "GPLv3"

# Different tokens needed for web access.
CHEMSPIDER_TOKEN = ""

#-------------------------------------------------------------------------------

class Compound(object):
  """ Class: Representing a chemical compund with its information details. """

  def __init__(self):
    """ Initial function for this class. """
    self._imageUrl = str()
    self._formula = str()
    self._smiles = str()
    self._name = str()
    self._csid = str()
    self._molWeight = str()
    self._database = str()

  @property
  def database(self):
    """ Holds the name of the database the compound was found. """
    return self._database

  @database.setter
  def database(self, database):
    self._database = database

  @database.deleter
  def database(self):
    del self._database

  @property
  def name(self):
    """ Holds the name of the compound. """
    return self._name

  @name.setter
  def name(self, name):
    self._name = name

  @name.deleter
  def name(self):
    del self._name

  @property
  def smiles(self):
    """ Holds the smiles information. """
    return self._smiles

  @smiles.setter
  def smiles(self, smiles):
    self._smiles = smiles

  @smiles.deleter
  def smiles(self):
    del self._smiles

  @property
  def csid(self):
    """ Holds the csid of the compund. """
    return self._csid

  @csid.setter
  def csid(self, csid):
    self._csid = csid

  @csid.deleter
  def csid(self):
    del self._csid

  @property
  def molWeight(self):
    return self._molWeight

  @molWeight.setter
  def molWeight(self, molWeight):
    self._molWeight = molWeight

  @molWeight.deleter
  def molWeight(self):
    del self._molWeight

  def __repr__(self):
    return "Compound(%s)" % (self._name)

#-------------------------------------------------------------------------------

class AbstractOnlineDatabase(metaclass = abc.ABCMeta):
  """ Class: Abstract class for database interface on the web. Every online
  database should have a class derived from this to be visible. """

  def __init__(self):
    """ Initial function for this class. """
    pass

  @abstractmethod
  def name(self):
    pass

  @abstractmethod
  def search(self, query):
    pass

#-------------------------------------------------------------------------------

class PubChem(AbstractOnlineDatabase):
  """ Class: Access to the web database of ChemPub. """

  def __init__(self):
    """ Initial function for this class. """
    super(AbstractOnlineDatabase, self)

  def name(self):
    """ Get the name of the online database. """
    return "PubChem"

  def search(self, query):
    """ Search compound and return it. """
    assert type(query) == str or type(query) == unicode, \
      "query not a string object"

    compounds = list()

    # Create the search url for access.
    resturl = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids" \
      "/XML" % query

    # Get the response from the url.
    response = urllib.request.urlopen(resturl)

    # Parse the root xml data from the response.
    root = ET.parse(response).getroot()

    # Return a list of cid elements.
    for child in root[::-1]:
      compounds.append(self.getCompoundFromCID(child.text))
      break

    return compounds

  def getCompoundFromCID(self, cid):
    """ Get a compund object with informations. """
    assert type(cid) == str or type(cid) == unicode, \
      "query not a string object"

    compound = Compound()
    compound.database = self.name()
    compound.csid = cid

    resturl = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/" \
      "record/XML" % cid

    # Get the response form the url.
    response = urllib.request.urlopen(resturl)

    # Parse the root xml data form the response.
    root = ET.parse(response).getroot()
    props = [element for element in root.iter() if element.tag == "{http://www.ncbi.nlm.nih.gov}PC-InfoData"]

    # Parse all informations about the compound.
    for child in props:
      label = None
      for item in child.iter():
        if item.tag == "{http://www.ncbi.nlm.nih.gov}PC-Urn_label":
          label = item
          break

      value = None
      for item in child.iter():
        if item.tag == "{http://www.ncbi.nlm.nih.gov}PC-InfoData_value_sval":
          value = item
          break

      if label.text == "IUPAC Name":
        compound.name = value.text

      if label.text == "SMILES":
        compound.smiles = value.text

    return compound

#-------------------------------------------------------------------------------

class ChemSpiderOD(AbstractOnlineDatabase):
  """ Class: Access to the web database of ChemSpider. """

  def __init__(self):
    """ Initial function for this class. """
    super(AbstractOnlineDatabase, self)

  def name(self):
    """ Get the name of the online database. """
    return "ChemSpider"

  def search(self, query):
    """ Search compound and return it. """
    assert type(query) == str or type(query) == unicode, \
      "query not a string object"

    compounds = list()

    # Create the search url for access.
    searchurl = "http://www.chemspider.com/Search.asmx/" \
      "SimpleSearch?query=%s&token=%s" % (urllib.parse.quote(query), \
       CHEMSPIDER_TOKEN)

    # Get the response form the url.
    response = urllib.request.urlopen(searchurl)

    # Parse the root xml data from the response.
    root = ET.parse(response).getroot()

    # Return a list of csid elements.
    for child in root:
      compounds.append(self.getCompoundFromCSID(child.text))

    return compounds

  def getCompoundFromCSID(self, csid):
    """ Get a compund object with informations. """
    assert type(csid) == str or type(csid) == unicode, \
      "query not a string object"

    compound = Compound()
    compound.database = self.name()

    apiurl = "http://www.chemspider.com/MassSpecAPI.asmx/" \
      "GetExtendedCompoundInfo?CSID=%s&token=%s" % (csid, CHEMSPIDER_TOKEN)

    # Get the response form the url.
    response = urllib.request.urlopen(apiurl)

    # Parse the root xml data form the response.
    root = ET.parse(response).getroot()

    # Parse all informations about the compound.
    for child in root:
      tag = child.tag
      text = child.text
      if tag == "{http://www.chemspider.com/}SMILES":
        compound.smiles = text
      elif tag == "{http://www.chemspider.com/}CSID":
        compound.csid = text
      elif tag == "{http://www.chemspider.com/}MolecularWeight":
        compound.molWeight = text
      elif tag == "{http://www.chemspider.com/}CommonName":
        compound.name = text

    return compound

#-------------------------------------------------------------------------------

class ChemInfo(object):
  """ Class: Main class for functionalitiy of pyChemInfo. """

  def __init__(self):
    """ Initial function for this class. """
    pass

  def search(self, query, first = False):
    """ Search for query and return all compounds from different online
    databases.

    query: The search string.
    first: If true, take the first element information found. """
    compounds = list()

    for odb in self.getOnlineDatabases():
      compounds += odb().search(query)

      # Continue if found a compound.
      if first is True and len(compounds) > 0:
        continue

    return compounds

  def getOnlineDatabases(self):
    return AbstractOnlineDatabase.__subclasses__()

#-------------------------------------------------------------------------------

def getODBClasses():
  AbstractOnlineDatabase.__subclasses__()

def main():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hfo:n:lv", ["help", "output=", \
        "names=", "list-databases"])
  except getopt.GetoptError as error:
    print(str(error))
    usage()
    sys.exit(2)

  output = None
  forceOutput = False
  names = None
  verbose = False
  for o, a in opts:
    if o == "-v":
      verbose = True
    elif o in ("-f"):
      forceOutput = True
    elif o in ("-l", "--list-databases"):
      print("List of available Online Database Modules: ", end="")
      print(",".join([chemInfo().name() for chemInfo in \
        ChemInfo().getOnlineDatabases()]))
      sys.exit(0)
    elif o in ("-n", "--names"):
      names = a.split(",")
    elif o in ("-h", "--help"):
      usage()
      sys.exit()
    elif o in ("-o", "--output"):
      output = a
    else:
      assert False, "unhandle option"

  if names is None:
    print("You are pleased to specify one or more name(s).")
    #usage()
    sys.exit(2)

  # Search for all names and compounds.
  compounds = None
  for name in names:
    chemInfo = ChemInfo()
    compounds = chemInfo.search(name, True)

    # Any compounds found?
    if len(compounds) < 1:
      print("Cannot find '%s' in %s." % (name, odb.name()))
      continue

  print("")
  for compound in compounds:
    # Print informations about the first compound.
    print("---------")
    print("Database: %s" % compound.database)
    print("Name:     %s" % compound.name)
    print("CSID:     %s" % compound.csid)
    print("SMILES:   %s" % compound.smiles)
  print("---------\n")

  # Write the output to a file if file name was given.
  if output is not None:
    if os.path.exists(output) is True and forceOutput is False:
      print("The file already exist. I won't write any data.")
    else:
      outFile = open(output, "w")
      for compound in compounds:
        # Print informations about the first compound.
        outFile.write("---------\n")
        outFile.write("Database: %s\n" % compound.database)
        outFile.write("Name:     %s\n" % compound.name)
        outFile.write("CSID:     %s\n" % compound.csid)
        outFile.write("SMILES:   %s\n" % compound.smiles)
      outFile.write("---------")
      outFile.close()

# Main entry point from command line calling.
if __name__ == "__main__":
  main()
