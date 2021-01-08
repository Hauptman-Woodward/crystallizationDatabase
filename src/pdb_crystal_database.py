import sys, json, pickle, operator, traceback, os, csv, itertools
from misc_functions import loadJson, writeJson, printList, fileToList, listToFile, getKey
from time import sleep
from collections import OrderedDict
from pathlib import Path
import xml.etree.ElementTree as etree
from xml.dom import minidom

# Make sure directories exist
if not os.path.exists("Output"):
    os.makedirs("Output")
if not os.path.exists("Structures"):
    os.makedirs("Structures")

# Create Path objects for directories
INPUT_DIR = Path("Input/")
STRUCTURE_DIR = Path("Structures/")
OUTPUT_DIR = Path("Output/")

# Input Files
COMPOUND_DICTIONARY_FILE = INPUT_DIR / "compound_dictionary.json"
SMILES_DICTIONARY_FILE = INPUT_DIR / "smiles_dictionary.json"

STOP_WORDS_FILE = INPUT_DIR / "stop_words.json"
UNKNOWN_LIST_FILE = INPUT_DIR / "unknown_list.json"
LOWERCASE_REPLACEMENT_FILE = INPUT_DIR / "replacementLowercase.json"
SENSITIVE_REPLACEMENT_FILE = INPUT_DIR / "replacementSensitive.json"
MIXTURES_FILE = INPUT_DIR / "mixture_compounds.json"
CLASSIFICATION_FILE = INPUT_DIR / "classification_dictionary.json"

# Structure Files - serialized binary files containing a list of Structure objects
STRUCTURES_FILE = STRUCTURE_DIR / "structures.pkl" # The database file. Must be placed in proper location
SENSIBLE_STRUCTURES_FILE = STRUCTURE_DIR / "sensible_structures.pkl" # Output only, not required

CSV_FILE = STRUCTURE_DIR / "sensible_structures.csv"
XML_FILE = STRUCTURE_DIR / "sensible_structures.xml"

# Output Files
DETAILS_FILE = OUTPUT_DIR / "details.txt"
SENSIBLE_DETAILS_FILE = OUTPUT_DIR / "sensible_details.txt"
NON_SENSIBLE_DETAILS_FILE = OUTPUT_DIR / "nonsensible_details.txt"
UNKNOWN_DETAILS_FILE = OUTPUT_DIR / "unknown_details.txt"

COMPOUND_FREQUENCY_FILE = OUTPUT_DIR / "compound_frequency.txt"
UNKNOWN_FREQUENCY_FILE = OUTPUT_DIR / "unknown_frequency.txt"
PENDING_FREQUENCY_FILE = OUTPUT_DIR / "pending_frequency.txt"
SET_FREQUENCY_FILE = OUTPUT_DIR / "set_frequency.txt"

COMPOUND_FREQUENCY_CSV_FILE = OUTPUT_DIR / "compound_frequency.csv"
UNKNOWN_FREQUENCY_CSV_FILE = OUTPUT_DIR / "unknown_frequency.csv"
PENDING_FREQUENCY_CSV_FILE = OUTPUT_DIR / "pending_frequency.csv"
SET_FREQUENCY_CSV_FILE = OUTPUT_DIR / "set_frequency.csv"


COMPRESSED_DICTIONARY_FILE = OUTPUT_DIR / "compressed_dictionary.json"

# Load lists and dictionaries from files
print("Loading input files...")
try:
    unknownList = loadJson(UNKNOWN_LIST_FILE)
    compoundDictionary = loadJson(COMPOUND_DICTIONARY_FILE)
    smilesDictionary = loadJson(SMILES_DICTIONARY_FILE)
    lowercaseReplacement = loadJson(LOWERCASE_REPLACEMENT_FILE, object_pairs_hook=OrderedDict)
    sensitiveReplacement = loadJson(SENSITIVE_REPLACEMENT_FILE, object_pairs_hook=OrderedDict)
    mixturesDictionary = loadJson(MIXTURES_FILE)
    classificationDictionary = loadJson(CLASSIFICATION_FILE)
    STOP_WORDS = set(loadJson(STOP_WORDS_FILE)) - {'m'} - {'am'}
except FileNotFoundError as notFound:
    print(notFound)
    print("ERROR: The file {} cannot be found. Verify that it is in the proper directory.".format(notFound.filename))

# Global variables set by exportOutputFiles()
compoundFrequency = {}
unknownFrequency = {}
pendingFrequency = {}
sensibleStructureList = []
nonsensibleStructureList = []

# Compounds that contain numbers (eg jeffamine 600)
NUMBERED_COMPOUNDS = ["jeffamine", "propoxylate", "polypropylene", "ndsb"]

class Structure:

    # Used for legacy version of removing protein solution
    # RESEVOIR_INDICATOR_WORDS = {"reservoir", "resevoir", "crystallization", "precipitant", "well", "solution", "cocktail"}

    def __init__(self, pdbid, pmcid, details, compounds, pH, temperature, method, sequences, resolution):
        self.pdbid = pdbid # String
        self.pmcid = pmcid # String
        self.details = details # crystallization details, String
        self.compounds = compounds # List of compounds followed by concentration
        self.pH = pH # float
        self.temperature = temperature # float, in K
        self.method = method # crystallization method, String
        self.sequences = sequences # List
        self.resolution = resolution # float, in angstroms
        self.compounds = []
        #self.parseDetails() # Sets the compounds list

    def __str__(self):
        return("""
        PDB: {0}
        PMCID: {1}
        Compounds: {2}
        pH: {3}
        Temperature: {4}
        Method: {5}
        Resolution: {6}
        """.format(self.pdbid, self.pmcid, self.compounds, self.pH, self.temperature, self.method, self.resolution))

    def parseDetails(self, detailsString=None, debug=False):
        """Parses the details string and returns a list of compounds, followed by concentration, or None if conc. is not found
        Format of compounds = ['compound name', '100', 'another compound name', '45%']
        By default (if detailsString=None), the input will be self.details, but
        detailsString can be set to another string instead
        If debug=True, then this function will print the words after every step
        Also sets the compounds list of the structure
        Returns None if failed or if details are unavailible
        """
        global nltk

        # Make sure NLTK is imported
        importNLTK()

        if (detailsString==None):
            details = self.details
        else:
            details = detailsString

        if details == None:
            return None

        if debug:
            print("Debugging PDB {}...\n".format(self.pdbid))
            print("Raw details:\n"+details+"\n")

        # Add spaces after commas before numbers (eg "sodium acetate,4% PEG4K" --> "sodium acetate, 4% PEG4K")
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(1, len(details)-1):
                if details[j] == "," and details[j-1].isalpha() and isNumber(details[j+1]):
                    # Insert a space
                    details = details[:j+1] + " " + details[j+1:]
                    runAgain = True
                    break

        if debug:
            print("Dealing with lack of space after commas:\n"+details+"\n")

        details = wordReplacement(details, sensitiveReplacement, debug=debug) # Replace with case sensitivity
        details = wordReplacement(details.lower(), lowercaseReplacement, debug=debug) # Replace without case sensitivity

        if debug:
            print("Word replacement:\n"+details+"\n")

        # Remove all text after 'cryo'
        cryoLocation = details.find("cryo")
        if cryoLocation != -1:
            details = details[:cryoLocation]

        # Removed all text after 'soaked'
        soakedLocation = details.find("soaked")
        if soakedLocation != -1:
            details = details[:soakedLocation]

        # Remove all text before "reservoir solution" to include only resevoir soltuion
        reservoirLocation = details.find("reservoir solution")
        if reservoirLocation != -1:
            details = details[reservoirLocation:]
            # Remove protein solution if it comes after reservoir solution
            proteinLocation = details.find("protein solution")
            if proteinLocation != -1:
                details = details[:proteinLocation]

        if debug:
            print("Isolate resevoir solution:\n"+details+"\n")

        words = nltk.word_tokenize(details)

        if debug:
            print("Tokenized words:\n"+str(words)+"\n")

        # Remove protein solution - legacy, obsolete
        # if len(words) > 1 and words[0] == "protein" and words[1] in Structure.RESEVOIR_INDICATOR_WORDS:
        #     reservoirSolutionFound = False
        #     for j in range(0, len(words)-1):
        #         if ((words[j] in Structure.RESEVOIR_INDICATOR_WORDS) and (words[j+1] in Structure.RESEVOIR_INDICATOR_WORDS)):
        #             words = words[j:]
        #             reservoirSolutionFound = True
        #             break
        #     if (not reservoirSolutionFound):
        #         return None # Only protein solution

        if debug:
            print("Remove protein solution:\n"+str(words)+"\n")

        # Fix commas between numbers (turn 6000,10 in to 6000 , 10)
        runAgain = True
        while (runAgain):
            runAgain = False
            for j in range(0, len(words)):
                if "," in words[j] and isNumber(words[j].replace("%","").replace("m","").replace("k","")):
                    numbers = words[j].split(",")
                    if len(numbers) == 2 and (len(numbers[1]) != 3 or len(numbers[0]) > 3 or "." in numbers[1]) and numbers[0] != "" and numbers[1] != "":
                        runAgain = True
                        words[j] = numbers[0]
                        words.insert(j+1, ",")
                        words.insert(j+2, numbers[1])
                        if len(numbers[1]) > 2 and numbers[1][-2:] == "mm":
                            words[j+2] = numbers[1][:-2]
                            words.insert(j+3, "mm")
                        elif numbers[1][-1:] == "m" or numbers[1][-1:] == "%":
                            words[j+2] = numbers[1][:-1]
                            words.insert(j+3, numbers[1][-1])
                        break

        if debug:
            print("Fix commas:\n"+str(words)+"\n")

        # Process micromolar (uM) concentrations, add to number like: "35 um" --> "35um"
        runAgain = True
        while (runAgain):
            runAgain = False
            for j in range(1, len(words)):
                if words[j] == "um" and isNumber(words[j-1]):
                    words[j-1] = words[j-1]+"um"
                    del words[j]
                    runAgain = True
                    break

        if debug:
            print("Process uM concentrations:\n"+str(words)+"\n")

        # Remove pH and temperature
        runAgain = True
        while (runAgain): # pH
            runAgain = False
            for j in range(0, len(words)):
                if j+1 < len(words) and words[j] == "ph" and isNumber(words[j+1]):
                    del words[j:j+2]
                    runAgain = True
                    break
                elif len(words[j]) > 2 and words[j][:2] == "ph" and isNumber(words[j][2:].replace("=","").replace(" ","")):
                    del words[j]
                    runAgain = True
                    break
                elif j+2 < len(words) and words[j] == "ph" and words[j+1] == "=" and isNumber(words[j+2]):
                    del words[j:j+3]
                    runAgain = True
                    break
        runAgain = True
        while (runAgain): # Temperature
            runAgain = False
            for j in range(1, len(words)):
                if (words[j][-1] == 'k' and len(words[j]) > 3 and isNumber(words[j][:-1])):
                    del words[j]
                    runAgain = True
                    break
                elif (words[j][-1] == 'k' and len(words[j]) > 1 and isNumber(words[j][:-1]) and words[j-1] == "temperature"):
                    del words[j]
                    runAgain = True
                    break
                elif words[j] == 'k' and isNumber(words[j-1]):
                    temp = float(words[j-1].replace(",",""))
                    if temp > 200 and temp < 400:
                        del words[j-1:j+1]
                        runAgain = True
                        break

        if debug:
            print("Remove pH and temperature:\n"+str(words)+"\n")

        # Process compounds with numbers in them
        for j in range(0, len(words)-1):
            if words[j] in NUMBERED_COMPOUNDS and isNumber(words[j+1]):
                words[j] = words[j] + " " + words[j+1]
                del words[j+1]
                break

        # Process PEG and MPEG compounds
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(0, len(words)):
                if (words[j][:3] == "peg" or words[j][:4] == "mpeg"):
                    pegNumber = ""
                    if j+1<len(words) and words[j+1] in STOP_WORDS:
                        runAgain = True
                        del words[j+1]
                        break
                    elif len(words[j]) > 3 and words[j][:3] == "peg": # If the number is already attached, eg PEG5000
                        if j+1<len(words) and words[j+1] == "mme": # Ex. PEG3000 MME
                            pegNumber = "MME " + words[j][3:].replace(",","").replace("-", " ")
                            del words[j+1]
                        else:
                            pegNumber = words[j][3:].replace("-", " ").replace(",","")
                    elif len(words[j]) > 4 and words[j][:4] == "mpeg": #MPEG4000
                        pegNumber = words[j][4:].replace("-", " ").replace(",","")
                    elif j+2 < len(words) and (words[j+1] == "mme" or words[j+1] == "monomethyl") and words[j+2] not in STOP_WORDS: # ex. PEG MME 5000
                        if isNumber(words[j+2].replace("k","000")):
                            pegNumber = "MME " + words[j+2].replace(",","")
                        else:
                            pegNumber = "MME"
                        del words[j+2]
                        del words[j+1]
                    elif j+2 < len(words) and words[j+2] == "mme" and words[j+1] not in STOP_WORDS: # ex. PEG 5000 MME
                        if isNumber(words[j+1].replace("k","000")):
                            pegNumber = "MME " + words[j+1].replace(",","")
                        else:
                            pegNumber = "MME"
                        del words[j+2]
                        del words[j+1]
                    elif j+1<len(words) and words[j+1] not in STOP_WORDS: # Ex PEG 5000
                        if isNumber(words[j+1].replace("k","000")):
                            pegNumber = words[j+1].replace(",","")
                            del words[j+1]
                    if len(pegNumber) > 3 and pegNumber[-3:] == "mme":
                        pegNumber = "MME " + pegNumber[:-3]
                    pegNumber = pegNumber.replace("000k","000").replace("k","000").replace("-","")
                    if words[j][:4] == "mpeg":
                        words[j] = "MPEG "+pegNumber
                    else:
                        words[j] = "PEG "+pegNumber
                    runAgain = True
                    break

        if debug:
            print("Process PEG and MPEG compounds:\n"+str(words)+"\n")

        # Remove terminal periods from numbers
        for j in range(0, len(words)):
            if isNumber(words[j]) and words[j][-1] == ".":
                words[j] = words[j][:-1]

        if debug:
            print("Remove terminal periods from numbers:\n"+str(words)+"\n")

        # Put together numbers separated by "-"
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(1, len(words)-1):
                if words[j] == "-" and isNumber(words[j-1]):
                    if isNumber(words[j+1]) or (words[j+1][-1] == "m" and isNumber(words[j+1][:-1])) or (len(words[j+1]) > 2 and words[j+1][-2] == "mm" and isNumber(words[j+1][:-2])):
                        words[j-1] = printList(words[j-1:j+2], "")
                        del words[j:j+2]
                        runAgain = True
                        break

        if debug:
            print("Combine numbers separated by '-':\n"+str(words)+"\n")

        # Average ranges of numbers
        for j in range(0, len(words)):
            if '-' in words[j]:
                if averageNumberString(words[j]) != None: # If there is no letter in front of the numbers
                    words[j] = str(averageNumberString(words[j]))
                elif averageNumberString(words[j][:-1]) != None: # If there is one letter in front of the numbers (e.g. 'm' or '%')
                    words.insert(j+1, words[j][-1])
                    words[j] = str(averageNumberString(words[j][:-1]))
                elif averageNumberString(words[j][:-2]) != None: # If there are two letters in front of the numbers (i.e. Mm)
                    words.insert(j+1, words[j][-2:])
                    words[j] = str(averageNumberString(words[j][:-2]))

        # Remove temperature again
        runAgain = True
        while (runAgain): # Temperature
            runAgain = False
            for j in range(1, len(words)):
                if (words[j][-1] == 'k' and len(words[j]) > 3 and isNumber(words[j][:-1])):
                    del words[j]
                    runAgain = True
                    break
                elif (words[j][-1] == 'k' and len(words[j]) > 1 and isNumber(words[j][:-1]) and words[j-1] == "temperature"):
                    del words[j]
                    runAgain = True
                    break
                elif words[j] == 'k' and isNumber(words[j-1]):
                    temp = float(words[j-1].replace(",",""))
                    if temp > 200 and temp < 400:
                        del words[j-1:j+1]
                        runAgain = True
                        break

        if debug:
            print("Average ranges of numbers and remove temperature again:\n"+str(words)+"\n")

        if len(words) == 0:
            return None

        # Process concentrations
        runAgain = True
        while(runAgain):
            for j in range(0, len(words)):
                if words[j][-1] == "m" and isNumber(words[j][:-1]): # Units are M
                    words[j] = str(float(words[j][:-1].replace("%","").replace(",",""))*1000)
                    runAgain = True
                    break
                elif words[j][-2:] == "mm" and isNumber(words[j][:-2]): # Units are mM
                    words[j] = words[j][:-2]
                    runAgain = True
                    break
                elif words[j][-2:] == "um" and isNumber(words[j][:-2]): # Units are uM
                    words[j] = str(float(words[j][:-2].replace("%","").replace(",",""))/1000)
                    runAgain = True
                    break
                elif isNumber(words[j]) and j+1 < len(words): # Conc. is separated from number by space, e.g. "0.1 M"
                    if words[j+1] == "m":
                        del words[j+1]
                        words[j] = str(float(words[j].replace(",",""))*1000)
                        runAgain = True
                        break
                    elif words[j+1] == "mm":
                        del words[j+1]
                        words[j] = str(float(words[j].replace(",","")))
                        runAgain = True
                        break
                    elif words[j+1] == "%":
                        if j+2 < len(words) and (words[j+2] == "w/v" or words[j+2] == "v/v"):
                            units = "% " + words[j+2]
                            del words[j+2]
                        else:
                            units = "%"
                        del words[j+1]
                        words[j] = str(words[j]+units)
                        runAgain = True
                        break
                else:
                    runAgain = False

        if debug:
            print("Process concentrations:\n"+str(words)+"\n")

        # Remove excess punctuation
        for j in range(0, len(words)):
            if words[j] in {",", "."}:
                if all(word in {",", "."} for word in words[j:]):
                    del words[j:]
                    break

        if debug:
            print("Remove excess punctuation:\n"+str(words)+"\n")

        # Combline parts of chemical compounds
        runAgain = True
        while runAgain:
            runAgain = False
            for j in range(0, len(words)-1):
                if isCompound(words[j]) and isCompound(words[j+1]) and words[j][:4] != "PEG " and words[j+1][:4] != "PEG " and words[j][:4] != "MPEG " and words[j+1][:4] != "MPEG ":
                    k = j
                    chemicalCompound = []
                    while(k < len(words) and isCompound(words[k])):
                        chemicalCompound.append(words[k])
                        k += 1
                    del words[j+1:k]
                    words[j] = printList(chemicalCompound, " ")
                    runAgain = True
                    break

        if debug:
            print("Combine parts of chemical compounds:\n"+str(words)+"\n")

        # Remove periods from compounds
        for j in range(0, len(words)):
            if isCompound(words[j]):
                words[j] = words[j].replace(".","")

        if debug:
            print("Remove periods from compounds:\n"+str(words)+"\n")

        # Combine w/v or v/v with non-percent nubers, and assume it's supposed to be percent
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(0, len(words)-1):
                if isNumber(words[j]) and (words[j+1] == "w/v" or words[j+1] == "v/v"):
                    words[j] = words[j] + "% " + words[j+1]
                    del words[j+1]
                    runAgain = True
                    break

        # If percent concentration is before and w/v is after
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(0, len(words)-2):
                if isPercent(words[j]) and isCompound(words[j+1]) and (words[j+2] == "w/v" or words[j+2] == "v/v"):
                    words[j] = words[j] + " " + words[j+2]
                    del words[j+2]
                    runAgain = True
                    break

        if debug:
            print("Parse w/v and v/v:\n"+str(words)+"\n")

        # Determine if concentration or compounds come first
        concentrationBeforeCompound = True
        for j in range(0, len(words)-1):
            if isConcentraton(words[j]) and isCompound(words[j+1]) and getKey(words[j+1]) in compoundDictionary:
                concentrationBeforeCompound = True
                break
            elif isConcentraton(words[j+1]) and isCompound(words[j]) and getKey(words[j]) in compoundDictionary:
                concentrationBeforeCompound = False
                break

        # Remove stop words
        words = [w for w in words if w not in STOP_WORDS]

        if debug:
            print("Remove stop words:\n"+str(words)+"\n")

        # Add compounds to list
        compounds = []
        if len(words) > 0 and words[0] != "ERROR":
            runAgain = True
            while(runAgain):
                runAgain = False
                for j in range(0, len(words)):
                    if isCompound(words[j]):
                        if concentrationBeforeCompound and j > 0 and isConcentraton(words[j-1]): # If a concentration is found (before compound)
                            compound = words[j]
                            concentration = words[j-1]
                            # If w/v or v/v comes after compound
                            if j+1<len(words) and "%" in concentration and (words[j+1] == "w/v" or words[j+1] == "v/v"):
                                concentration = concentration + " " + words[j+1]
                                del words[j-1:j+2]
                            else:
                                del words[j-1:j+1]
                        elif not concentrationBeforeCompound and j+1 < len(words) and isConcentraton(words[j+1]): # if a concentration is found (after compound)
                            compound = words[j]
                            concentration = words[j+1]
                            del words[j:j+2]
                        else: # No concentration found
                            compound = words[j]
                            concentration = None
                            del words[j]
                        compounds.append(compound)
                        compounds.append(concentration)
                        runAgain = True
                        break
        self.compounds = compounds

        if debug:
            print("Extract compounds:\n"+str(compounds)+"\n")

        return compounds

    def standardizeNames(self):
        """Standardizes the name of an individual structure"""
        runAgain = True
        while(runAgain):
            runAgain = False
            for i in range(0, len(self.compounds), 2):
                if getKey(self.compounds[i]) in compoundDictionary:
                    self.compounds[i] = compoundDictionary[getKey(self.compounds[i])]

                    # Parse multiple compounds delimited with " / "
                    if " / " in self.compounds[i]:
                        multipleCompounds = self.compounds[i].split(" / ")
                        concentration = self.compounds[i+1]
                        for c in multipleCompounds:
                            if c in set(compoundDictionary.values()):
                                self.compounds.append(c)
                                self.compounds.append(concentration)
                            else:
                                print("ERROR: Compound {} was found in multiple compound {} but is not in dictionary.".format(c, self.compounds[i]))
                                input = input("Add the key {} : {} to the dictionary? (y/n): ".format(getKey(c), c))
                                if input == "y":
                                    compoundDictionary[getKey(c)] = c
                                    print("Added key to dictionary")
                                    self.compounds.append(c)
                                    self.compounds.append(concentration)
                        del self.compounds[i:i+2]
                        runAgain = True
                        break

                    # Parse mixture compounds, that is,
                    # names of mixtures which refer to multiple chemical compounds
                    if self.compounds[i] in mixturesDictionary:
                        mixtureDict = mixturesDictionary[self.compounds[i]] # Dictionary of the compounds in the mixture
                        concentration = self.compounds[i+1]
                        if concentration != None:
                            findPercent = concentration.find("%")
                            percentUnits = ""
                            if findPercent != -1:
                                percentUnits = concentration[findPercent:] # Everything after and including '%'
                                concentration = concentration[:findPercent] # Everything berfore "%"
                        for c in mixtureDict:
                            if c in set(compoundDictionary.values()):
                                self.compounds.append(c)
                                if concentration != None:
                                    if percentUnits == "":
                                        newConcentration = str(float(concentration)*mixtureDict[c])
                                    else:
                                        # If concentration is a percent, then calculate millimolar concentration
                                        newConcentration = str(float(concentration)*mixtureDict[c]*10)
                                else:
                                    newConcentration = None
                                self.compounds.append(newConcentration)
                            else:
                                print("ERROR: Compound {} was found in mixture compound {} but is not in dictionary.".format(c, self.compounds[i]))
                                input = input("Add the key {} : {} to the dictionary? (y/n): ".format(getKey(c), c))
                                if input == "y":
                                    compoundDictionary[getKey(c)] = c
                                    print("Added key to dictionary")
                                    self.compounds.append(c)
                                    if concentration != None:
                                        if percentUnits == "":
                                            newConcentration = str(float(concentration)*mixtureDict[c])
                                        else:
                                            # If concentration is a percent, then calculate millimolar concentration
                                            newConcentration = str(float(concentration)*mixtureDict[c]*10)
                                    else:
                                        newConcentration = None
                                    self.compounds.append(newConcentration)
                        del self.compounds[i:i+2]
                        runAgain = True
                        break

                    # Remove duplicate compounds
                    if self.compounds[i] in self.compounds[:i:2]: # If compound has appeared before
                        locationOfFirstCompound = self.compounds.index(self.compounds[i], 0, i)
                        if self.compounds[i+1] == None and self.compounds[locationOfFirstCompound+1] != None: # If only the first compound has a concentrations
                            del self.compounds[i:i+2] # Remove the second compound
                            runAgain = True
                            break
                        else:
                            del self.compounds[locationOfFirstCompound:locationOfFirstCompound+2] # Remove first compound
                            runAgain = True
                            break

    def getXml(self):
        """Returns an elementTree tree root of the structure object"""
        root = etree.Element("structure")

        etree.SubElement(root, "pdbid").text = self.pdbid
        etree.SubElement(root, "pmcid").text = str(self.pmcid)
        etree.SubElement(root, "details").text = str(self.details)
        compoundsTag = etree.SubElement(root, "compounds")
        for i in range(0, len(self.compounds), 2):
            compoundTag = etree.SubElement(compoundsTag, "compound")
            etree.SubElement(compoundTag, "name").text = str(self.compounds[i])
            etree.SubElement(compoundTag, "concentration").text = str(self.compounds[i+1])
        etree.SubElement(root, "pH").text = str(self.pH)
        etree.SubElement(root, "temperature").text = str(self.temperature)
        etree.SubElement(root, "method").text = str(self.method)
        etree.SubElement(root, "resolution").text = str(self.resolution)
        sequencesTag = etree.SubElement(root, "sequences")
        for s in self.sequences:
            etree.SubElement(sequencesTag, "sequence").text = s

        return root

    def hasUnknown(self): # boolean
        """returns True if the structure has a compound in the unknownList"""
        for compound in self.compounds[::2]:
            if compound in unknownList:
                return True

    def printError(self, errorMessage, errorObject):
        """Prints an error message regarding a certain structure.
        errorMessage is a string which describes when the error occrued
        errorObject is the exception object which is caught"""
        print("--------------------\nERROR: {} for PDB ID {}.\n".format(errorMessage, self.pdbid))
        traceback.print_tb(errorObject.__traceback__)
        print(str(errorObject)+"\n")
        print("Crystallization Details: {}\n\nCompounds: {}\n--------------------\n".format(self.details, self.compounds))

def parseAllDetails(structureList, structureFile=None, searchString=None):
    """Reparses all of the details for a list of structures
    Should be called when the parseDetails function has been modified
    If search string is not None, then it will only parse Structures
    which have the search string in their details, making it faster
    If structureFile is specified, the structure list is saved to that file
    """
    global nltk
    # Make sure NLTK is imported
    if 'nltk' not in sys.modules:
        importNLTK()

    print("Parsing details with search string: \"{}\"...".format(searchString))
    count = 1

    if searchString != None:
        structureList = [structure for structure in structureList if searchString.lower() in structure.details.lower()]
        if structureList == []:
            print("No structures found with search string '{}'".format(searchString))

    print("Parsing details of {} structures...".format(len(structureList)))
    for structure in structureList:
        if count % 10000 == 0:
            print("Parsing structure {} of {}...".format(count, len(structureList)))
        try:
            structure.parseDetails()
        except Exception as e:
            structure.printError("Unable to parse details", e)
        count += 1

    if structureFile != None:
        print("Writing to structure file {}...".format(structureFile))
        writeStructures(structureList, structureFile)

def standardizeAllNames(structureList, structureFile=None): # void
    """Standardizes names of a list of compounds based on the compound Dictionary
    Also parses dictionary values which represent multiple compounds (e.g. acetic acid / sodium acetate)
    Also parses mixture compounds, which are mixtures of multiple compounds (e.g Molecular Dimensions Buffers)
    If structureFile is specified, the structure list is saved to that file
    """
    updateMiscDictionaries() # Also calls updateDictionary
    print("Standardizing chemical names...")
    for structure in structureList:
        try:
            structure.standardizeNames()
        except Exception as e:
            structure.printError("Unable to standardize compound names", e)
    if structureFile != None:
        print("Writing to structure file {}...".format(structureFile))
        writeStructures(structureList, structureFile)

def importNLTK():
    """Handles importing the NLTK Module"""
    global nltk
    if 'nltk' not in sys.modules:
        try:
            import nltk
            downloadPunkt()
        except ModuleNotFoundError:
            # Print traceback
            for line in traceback.format_stack():
                print(line.strip())

            print("ERROR: NLTK module not found. You must install the NLTK module in order to parse the details of structures. "
            "You may search through the database without NLTK, but you can not parse the details. "
            "For more information, see https://github.com/maxdudek/crystallizationDatabase/wiki/Using-the-Database-Script#Installing_NLTK")
            sys.exit()
    else:
        downloadPunkt()

def downloadPunkt():
    global nltk
    """Attempt to download the 'punkt' package from nltk"""
    try:
        import nltk
        nltk.data.find('tokenizers/punkt')
    except LookupError:
        try:
            print("Attempting to download 'punkt' package from nltk...")
            nltk.download('punkt')
            print("Successfully downloaded 'punkt' package.")
        except:
            # Print traceback
            for line in traceback.format_stack():
                print(line.strip())

            print("ERROR: The NLTK package 'punkt' could not be downloaded, and so you will be unable to parse details. For more information, "
            "see https://github.com/maxdudek/crystallizationDatabase/wiki/Using-the-Database-Script#Installing_NLTK")
            sys.exit()
    except ModuleNotFoundError:
        # Print traceback
        for line in traceback.format_stack():
            print(line.strip())

        print("ERROR: You should never get this error. If you do, it means that "
        "the script attempted to download the nltk 'punkt' package when nltk was not installed")
        sys.exit()

def updateMiscDictionaries(): # void
    """Adds newly recognized compounds in the compound dictionary as blank keys to the SMILES dictionary and to the classificationDictionary
    Also removes unessicary entries from these dictionaries"""
    updateDictionary()
    print("Updating misc. dictionaries...")

    # Add new compounds to dictionaries
    for compound in set(compoundDictionary.values()):
        if compound not in smilesDictionary:
            smilesDictionary[compound] = ""

    for compound in set(compoundDictionary.values()):
        if compound not in classificationDictionary:
            classificationDictionary[compound] = []

    # Remove mixture compounds from dictionaries
    for compound in set(compoundDictionary.values()):
        if " / " in compound or compound in mixturesDictionary:
            smilesDictionary.pop(compound, None)
            classificationDictionary.pop(compound, None)

    # Remove elements that are not in the compound dictionary
    for compound in set(smilesDictionary.keys()) | set(classificationDictionary.keys()):
        if compound not in compoundDictionary.values() and (compound in smilesDictionary or compound in classificationDictionary):
            smilesDictionary.pop(compound, None)
            classificationDictionary.pop(compound, None)

    # Autoclassify some compounds
    for compound in classificationDictionary:
        if classificationDictionary[compound] == []:
            if (compound[:4] == "PEG " or compound[:5] == "MPEG "):
                classificationDictionary[compound] = ['Polymer']
            if "(II)" in compound or "(III)" in compound or "(I)" in compound:
                classificationDictionary[compound] = ['Salt']

    # Save dictionaries to files
    writeJson(smilesDictionary, SMILES_DICTIONARY_FILE, indent=2, sort_keys=True)
    writeJson(classificationDictionary, CLASSIFICATION_FILE, indent=2, sort_keys=True)

def getSensibleStructures(structureList):
    """Returns a list of structures that make sense
    That is, all of their compounds are found in the dictionary"""
    print("Searching for sensible structures...")
    sensibleStructures = []
    for structure in structureList:
        isSensible = True
        if structure.compounds == []:
            isSensible = False
        for compound in structure.compounds[::2]:
            if compound not in set(compoundDictionary.values()):
                isSensible = False
                break
        if isSensible:
            sensibleStructures.append(structure)
    return sensibleStructures

def exportSensibleStructures(structureList, structureFile=SENSIBLE_STRUCTURES_FILE, detailsFile=SENSIBLE_DETAILS_FILE):
    """Exports a list of sensible structures to a serialized pickle file
    Also returns the list of sensible structures from getSensibleStructures
    Also exports the details of the sensible structures"""
    sensibleStructures = getSensibleStructures(structureList)
    writeStructures(sensibleStructures, structureFile)
    print("Retrieved {} sensible structures".format(len(sensibleStructures)))

    exportDetails(sensibleStructures, detailsFile)
    return sensibleStructures

def getUnknowns(structureList):
    """Takes a structure list, and returns a list of all structures which have an unknown compound"""
    print("Finding unknown compounds...")
    return [structure for structure in structureList if structure.hasUnknown()]

def exportDetails(structureList, outputFilename):
    """Exports the details of a list of structures to a text file"""
    print("Exporting details to \"{}\"...".format(outputFilename))
    outputList = []
    for structure in structureList:
        if structure.compounds != []:
            outputList.append(structure.pdbid + ": " + structure.details)
            outputList.append(structure.compounds)
    listToFile(outputList, outputFilename)

def getCompoundFrequencies(structureList, textFilename=None, csvFilename=None, mode="recognized"):
    """Takes a list of structures and returns a dictionary mapping compounds to their frequency
    If textFilename != None, it also outputs compound frequency to an easily readable txt file
    If csvFilename != None, then it also exports compound frequency to a TAB-DELIMITED csv file
    Mode = "recognized": only export compounds found in the dictionary
    Mode = "unknown": only export compounds found in unknownList
    Mode = "pending": only export compounds in neither dictionary nor unknown list, pending classification
    """
    global compoundDictionary
    global smilesDictionary

    print("Exporting compound frequencies to \"{}\"...".format(textFilename))
    outputDictionary = {"Total Compounds": 0}
    if mode == "recognized":
        for structure in structureList:
            for compound in structure.compounds[::2]:
                c = compound
                if getKey(c) in compoundDictionary:
                    c = compoundDictionary[getKey(c)]
                if c in compoundDictionary.values():
                    if c not in outputDictionary:
                        outputDictionary[c] = 1
                    else:
                        outputDictionary[c] += 1
                    outputDictionary["Total Compounds"] += 1

    if mode == "unknown":
        for structure in structureList:
            for compound in structure.compounds[::2]:
                if getKey(compound) in unknownList and getKey(compound) not in compoundDictionary:
                    if compound not in outputDictionary:
                        outputDictionary[compound] = 1
                    else:
                        outputDictionary[compound] += 1
                    outputDictionary["Total Compounds"] += 1

    if mode == "pending":
        for structure in structureList:
            for compound in structure.compounds[::2]:
                if getKey(compound) not in unknownList and getKey(compound) not in compoundDictionary and compound not in compoundDictionary.values():
                    if compound not in outputDictionary:
                        outputDictionary[compound] = 1
                    else:
                        outputDictionary[compound] += 1
                    outputDictionary["Total Compounds"] += 1

    # Export
    if textFilename != None:
        outputList = []
        for key, value in sorted(outputDictionary.items(), key=operator.itemgetter(1), reverse=True):
            outputList.append("{:30s}: {:20d}".format(key, value))
        listToFile(outputList, textFilename)

    if csvFilename != None:
        outputList = []
        for key, value in sorted(outputDictionary.items(), key=operator.itemgetter(1), reverse=True):
            outputList.append([key, value])
        with open(csvFilename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter="\t")
            writer.writerows(outputList)

    return outputDictionary

def getSetFrequencies(structureList, textFilename=None, csvFilename=None, requiredCompounds=None, subsetLength=None):
    """Takes a compound list and outputs a dictionary of the frequencies of each set of compounds
        The dictionary maps a frozenset of compounds to its frequency as an integer.
    If textFilename != None, then it will also export the dictionary to an easily readable text file
    If csvFilename != None, then it also exports compound frequency to a TAB-DELIMITED csv file
    If requiredCompounds != None, then it will only find sets that contain a specific compound
    If subsetLength != None, then teh function will calculate the frequency of SUBSETS of length 'subsetLength'
        For example, if subsetLength = 3, the function will calculate the frequency of the most common 3 element subsets
        This means that a frequency count for the set 'HEPES', 'ammonium sulfate', and 'PEG 3350' will include
        the frequency of ANY set which contains those three elements"""
    if textFilename != None or csvFilename != None:
        print("Exporting set frequencies to {} and {}...".format(textFilename, csvFilename))

    outputDictionary = {}

    if type(requiredCompounds) == str:
        structureList = [s for s in structureList if requiredCompounds in s.compounds[::2]]
    elif type(requiredCompounds) == list:
        for requiredCompound in requiredCompounds:
            structureList = [s for s in structureList if requiredCompound in s.compounds[::2]]

    for structure in structureList:
        compoundSet = frozenset(structure.compounds[::2])
        if compoundSet in outputDictionary:
            outputDictionary[compoundSet] += 1
        else:
            outputDictionary[compoundSet] = 1

    if subsetLength != None:
        newOutputDictionary = {}
        setOfAllSubsets = frozenset([])
        for s in outputDictionary: # Build a set of all subsets
            setOfAllSubsets = setOfAllSubsets | frozenset(itertools.combinations(s, subsetLength))

        # Clean up the set
        setOfAllSubsets = frozenset([frozenset(s) for s in setOfAllSubsets])

        count = 0
        for subset in setOfAllSubsets:
            if subset not in newOutputDictionary:
                count += 1
                if count % 5000 == 0:
                    print("Working on subset {} of {}: {}".format(count, len(setOfAllSubsets), str(subset)))
                listOfSupersets = [s for s in outputDictionary.keys() if s.issuperset(subset)]
                frequency = 0
                for superset in listOfSupersets:
                    frequency += outputDictionary[superset]
                newOutputDictionary[subset] = frequency
        outputDictionary = newOutputDictionary

    # Export
    if textFilename != None:
        outputList = []
        for key, value in sorted(outputDictionary.items(), key=operator.itemgetter(1), reverse=True):
            outputList.append("{:80s}: {:20d}".format(printList(list(key)), value))
        listToFile(outputList, textFilename)

    if csvFilename != None:
        outputList = []
        for key, value in sorted(outputDictionary.items(), key=operator.itemgetter(1), reverse=True):
            outputList.append([printList(list(key)), value])
        with open(csvFilename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter="\t")
            writer.writerows(outputList)

    return outputDictionary

def exportOutputFiles(structureList):
    """Just a simple way to export all of the output files
    Also sets a bunch of useful global variables"""
    global compoundFrequency, unknownFrequency, pendingFrequency, sensibleStructureList, nonsensibleStructureList
    sensibleStructureList = exportSensibleStructures(structureList)
    nonsensibleStructureList = list(set(structureList) - set(sensibleStructureList))
    exportDetails(structureList, DETAILS_FILE)
    exportDetails(getUnknowns(structureList), UNKNOWN_DETAILS_FILE)
    exportDetails(nonsensibleStructureList, NON_SENSIBLE_DETAILS_FILE)
    compoundFrequency = getCompoundFrequencies(sensibleStructureList, textFilename=COMPOUND_FREQUENCY_FILE, csvFilename=COMPOUND_FREQUENCY_CSV_FILE)
    unknownFrequency = getCompoundFrequencies(structureList, textFilename=UNKNOWN_FREQUENCY_FILE, csvFilename=UNKNOWN_FREQUENCY_CSV_FILE, mode="unknown")
    pendingFrequency = getCompoundFrequencies(structureList, textFilename=PENDING_FREQUENCY_FILE, csvFilename=PENDING_FREQUENCY_CSV_FILE, mode="pending")
    setFrequency = getSetFrequencies(sensibleStructureList, textFilename=SET_FREQUENCY_FILE, csvFilename=SET_FREQUENCY_CSV_FILE)
    exportXml(sensibleStructureList, XML_FILE)
    exportCsv(sensibleStructureList, CSV_FILE)

def exportCsv(structureList, outputFilename):
    """Exports a csv file of all of the information in a list of structures"""
    print("Exporting csv file to {}...".format(outputFilename))
    outputList = []
    for structure in structureList:
        row = [str(structure.pdbid), str(structure.pmcid), str(structure.method).replace(',',''), str(structure.resolution), str(structure.temperature), str(structure.pH)]
        for compound in structure.compounds:
            if compound != None:
                row.append(compound)
            else:
                row.append('None')
        while len(row) < 30: # Pad empty compound elements (6 columns + (12 compounds*2) = 30)
            row.append("---")
        for sequence in structure.sequences:
            row.append(sequence)
        outputList.append(row)

    try:
        with open(outputFilename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter="\t")
            writer.writerows(outputList)
    except PermissionError:
        print("Could not write to {}. Maybe the file is open in another program such as Excel?".format(outputFilename))

def exportXml(structureList, outputFilename):
    """Exports a list of structures to an xml file"""
    print("Exporting xml structure list to {}...".format(outputFilename))
    root = etree.Element("structures")
    for s in structureList:
        root.append(s.getXml())

    xmlString = minidom.parseString(etree.tostring(root)).toprettyxml(indent="   ")
    xmlLineList = xmlString.split("\n")

    try:
        with open(outputFilename, "w") as f:
            for line in xmlLineList:
                f.write(line)
                f.write("\n")
    except MemoryError:
        print("ERROR: Not enough memory availible to export xml to {}".format(outputFilename))

def getDatabaseSubset(structureList, pdbidList, sensibleOnly=True, structureFile=None):
    """Takes a list of PDB IDs and returns a list of structure files associated with them
    If sensibleOnly is True, then only 'sensible' structures will be included
    structureFile is an optional filename argument, if it is specified, then it will
    export the subset of structure objects to the structure file specified"""
    print("\nFetching database subset to export to {}...".format(structureFile))
    structureSubsetList = []
    pdbidList = list(set(pdbidList))
    for pdbid in pdbidList:
        structure = getStructure(structureList, pdbid)
        if structure != None:
            structureSubsetList.append(structure)

    print("Fetched {} structures for databse subset".format(len(structureSubsetList)))

    if sensibleOnly:
        structureSubsetList = getSensibleStructures(structureSubsetList)
        print("Fetched {} sensible structures from database subset".format(len(structureSubsetList)))

    if structureFile != None:
        print("Exporting structure subset to {}".format(structureFile))
        writeStructures(structureSubsetList, structureFile)

    return structureSubsetList

def updateDictionary():
    """Makes sure all values in the compound dictionary are also keys
    This prevents inconsistency in other functions which assume that this is the case"""
    for value in set(compoundDictionary.values()):
        if getKey(value) not in compoundDictionary:
            compoundDictionary[getKey(value)] = value
    writeJson(compoundDictionary, COMPOUND_DICTIONARY_FILE, indent=2)

def getStructure(structureList, pdbid): # Structure
    """Returns a specific structure in a list based on its pdbid"""
    for structure in structureList:
        if pdbid == structure.pdbid:
            return structure
    return None

def isNumber(s): # Boolean
    try:
        float(s.replace(",",""))
        return True
    except ValueError:
        return False # Boolean

def isConcentraton(s):
    return isNumber(s) or isPercent(s) # Boolean

def isPercent(s): # Boolean
    if (s[-3:] == "w/v" or s[-3:] == "v/v") and "%" in s:
        return True
    elif len(s) > 1 and s[-1] == "%":
        try:
            float(s[:-1])
            return True
        except ValueError:
            return False
    else:
        return False # Boolean

def isCompound(word): # boolean
    """Takes a string and returns true if the string is part of a chemical compound
    Returns false if the string is a comma or a number"""
    isPunctuation = word in [",", ".", ":", ";", "-"]
    return (not isConcentraton(word) and not isPunctuation and not word == "ERROR" and word not in STOP_WORDS)

def averageNumberString(numberRange): # float
    """Takes a string such as 20-25 and averages the two (or more) numbers (ex. "20-25" --> "22.5")
    Returns None if the split doesn't produce all numbers (ex. 20%-25%)
    """
    stringSplit = numberRange.split('-')
    if all(isNumber(number) for number in stringSplit):
        numbers = [float(i.replace(",","")) for i in stringSplit]
        return sum(numbers) / len(numbers)
    else:
        return None

def wordReplacement(s, replacementDictionary, debug=False):
    """Replaces certain strings to make parsing easier
    replacementDictionary: an OrderedDict dictionary which controls the replacements to make
    If debug is True, then every step is printed"""
    details = s
    ignoredSequence = "$*"
    for key, value in replacementDictionary.items():
        details = details.replace(key.replace(ignoredSequence, ""), value)
        if debug:
            print("Replacing '{}' with '{}':".format(key.replace(ignoredSequence, ""),value))
            print(details+"\n")
    return details

def loadStructures(structureFile=STRUCTURES_FILE): # list
    """Returns a list of structures from the pickled structure file"""
    print("Loading structures from file {}...".format(structureFile))
    with open(structureFile, "rb") as f:
        return pickle.load(open(structureFile, "rb"))

def writeStructures(structureList, structureFile, count=0):
    """Writes a list of structures to a pickle file
    count keeps track of how many times the function has had to wait to write"""
    if count > 5:
        print("ERROR: Permission denied {} times when trying to write structures to {}".format(count-1, structureFile))
        return None
    try:
        with open(structureFile, "wb") as f:
            pickle.dump(structureList, f)
            if count > 0:
                print("Successfully wrote structures")
            return True
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Writing to file...")
        with open(structureFile, "wb") as f:
            pickle.dump(structureList, f)
        sys.exit()
    except PermissionError:
        if count == 0:
            print("Permission denied when attempting to write structures. Waiting and trying again...")
        sleep(3)
        nextTry = writeStructures(structureList, structureFile, count=count+1)
        return nextTry

if __name__ == "__main__":
    # Insert your commands here
    structureList = loadStructures(STRUCTURES_FILE) # Must have a Structure File availible (see wiki)
    parseAllDetails(structureList, searchString=None, structureFile=STRUCTURES_FILE)
    standardizeAllNames(structureList, structureFile=STRUCTURES_FILE)
    exportOutputFiles(structureList)

    # Get most common partners
    # listOfMostFrequentCompounds = sorted(compoundFrequency.items(), key=operator.itemgetter(1), reverse=True)[1:11]
    # for c, freq in listOfMostFrequentCompounds:
    #     print("Getting '{}'".format(c))
    #     setFrequency = getSetFrequencies(sensibleStructureList, requiredCompounds=c, subsetLength=2)
    #     mostCommonSubsets = sorted(setFrequency.items(), key=operator.itemgetter(1), reverse=True)[:3]
    #     for p, f in mostCommonSubsets:
    #         print("\t{} - Frequency: {}".format(p, f))

    # Get max number of compounds
    # maxCompounds = 3
    # for structure in sensibleStructureList:
    #     if len(structure.compounds) / 2 > maxCompounds:
    #         maxStructure = structure
    #         maxCompounds = len(structure.compounds) / 2
    # print(maxCompounds)
    # print(maxStructure.details)
    # print(maxStructure.compounds)

    # Get number of PEG compounds
    # count = 0
    # for structure in sensibleStructureList:
    #     for c in structure.compounds[::2]:
    #         if c[:4] == "PEG ":
    #             count += 1
    # print("Number of PEG compounds: {}".format(count))

    # Debug parseDetails
    # getStructure(structureList, "3NNQ").parseDetails(debug=True)
