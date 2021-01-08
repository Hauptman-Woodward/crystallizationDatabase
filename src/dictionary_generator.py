import json, csv, sys, pickle, collections, os
from misc_functions import loadJson, writeJson, printList, getKey, listToFile, fileToList
from pdb_crystal_database import loadStructures, parseAllDetails, writeStructures, Structure
from pathlib import Path

# Make sure directories exist
if not os.path.exists("Output"):
    os.makedirs("Output")
if not os.path.exists("Structures"):
    os.makedirs("Structures")

# Create Path objects for directories
INPUT_DIR = Path("Input/")
STRUCTURE_DIR = Path("Structures/")
OUTPUT_DIR = Path("Output/")

# Configure input locations
COMPOUND_DICTIONARY_FILE = INPUT_DIR / "compound_dictionary.json"
UNKNOWN_LIST_FILE = INPUT_DIR / "unknown_list.json"
STOP_WORDS_FILE = INPUT_DIR / "stop_words.json"

STRUCTURES_FILE = STRUCTURE_DIR / "structures.pkl"

COMPRESSED_DICTIONARY_FILE = OUTPUT_DIR / "compressed_dictionary.json"

# Define input options
INPUT_SAME = "=" # Add the compound to the dictionary exactly as it appears (e.g. "sodium chloride" --> "sodium chloride")
INPUT_UNKNOWN = "unknown" # Add the compound to the unknownList
INPUT_STOP_WORDS = "sw" # Add ALL WORDS in the compound to the list of stop words
# ex. if the compound is "well plate", both "well" and "plate" are added to the list
INPUT_ADD_STOP_WORD = "add" # Adds a single stop word to the stopWord list
# ex. "add well" will append "well" to the stopWords list
INPUT_PASS = "pass" # Skips the current compound, does nothing to it
INPUT_UNDO = "u" # Undo - NOTE: undo is unreliable and should only be used for single dictionary changes
# To undo a mistake safely, EXIT THE SCRIPT and make the change in the text/json files directly
INPUT_SAVE = "save" # Save files - currently done automatically after every change
INPUT_QUIT = "quit" # Exit the script
INPUT_QUIT_WITHOUT_SAVING = "quit no save" # Exit the script without saving files

# Load input files
print("Loading input files for dictionary generator...")
try:
    compoundDictionary = loadJson(COMPOUND_DICTIONARY_FILE)
except FileNotFoundError:
    print("The compound dictionary file specified ({}) was not found. A blank dictionary file will be created at the specified location".format(COMPOUND_DICTIONARY_FILE))
    compoundDictionary = {}

try:
    stopWords = loadJson(STOP_WORDS_FILE)
except FileNotFoundError:
    print("The stop words file specified ({}) was not found. A blank file will be created at the specified location".format(STOP_WORDS_FILE))
    stopWords = []

try:
    unknownList = loadJson(UNKNOWN_LIST_FILE)
except FileNotFoundError:
    print("The unknown list file specified ({}) was not found. A blank file will be created at the specified location".format(UNKNOWN_LIST_FILE))
    unknownList = []

passList = [] # A temporary list to allow the user to temporarily skip a compound

def getCompoundList(structureList, sortedByFrequency=True, useGetKey=False): # list
    """Takes in a list of Structure Objects and returns a list of just the compound names
    If sortedByFrequency is True, the compound names are sorted by frequency
    If useGetKey is true, then the compounds are turned into dictionary keys using getKey(), and then sorted (eg "na-cl" --> "nacl")
    """
    print("Getting compound list...")
    compoundList = []
    for structure in structureList:
        compoundList.extend(structure.compounds[0::2])
    if useGetKey:
        compoundList = [getKey(compound) for compound in compoundList if compoundList]
    if sortedByFrequency:
        counts = collections.Counter(compoundList)
        compoundList = sorted(compoundList, key=lambda x: -counts[x])
    return compoundList

def getCompressedDictionary(dictionary, filename=None):
    """Takes a dictionary and returns a 'compressed' dictionary
    where every unique value in the original dictionary maps to a list of all keys which map to that value
    If filename is not None, then the output will be written to the specified json file"""
    outputDictionary = {}
    for value in set(dictionary.values()):
        outputDictionary[value] = []

    for key in dictionary:
        outputDictionary[dictionary[key]].append(key)

    if filename != None:
        writeJson(outputDictionary, filename, indent=2, sort_keys=True)

    return outputDictionary

def printRecognizedCompounds(compoundList):
    """Prints out how many compounds are recognized, out of the total compounds"""
    count = 0
    for compound in compoundList:
        if compound in compoundDictionary:
            count += 1
    print("{} out of {} compounds recognized".format(count, len(compoundList)))

def removeStopWords(s, stopWords):
    """Takes a string and returns a string with stop words removed
    stopWords = A list of words to remove"""
    words = s.split(" ")
    words = [word for word in words if word not in stopWords]
    return printList(words, " ")

def saveFiles():
    """Saves all of the lists and dictionaries to their respective files"""
    writeJson(compoundDictionary, COMPOUND_DICTIONARY_FILE, indent=2)
    writeJson(unknownList, UNKNOWN_LIST_FILE, indent=2)
    writeJson(stopWords, STOP_WORDS_FILE, indent=2)
    print("Files saved")

def generateDictionary(compoundList, autoSave=True, autoAdd=True): # dictionary
    """ Iterates through a list and substitutes elements based on a dictionary
    If no key is found for the element, the user is prompted to enter an entry
    See INPUT definitions above for more options
    If autoSave is True, then the files will save after every step
    If autoAdd is True, the dictionary will automatically add an identical key for every value entered
        For example, if "nacl" is mapped to "sodium chloride", the key "sodiumchloride" is also added to map to "sodium chloride"
    """
    print("Beginning dictionary generation (may take a minute)...")
    history = [] # A list of modified indeces, in order to UNDO
    i = 0
    while(i < len(compoundList)+1):
        if i < len(compoundList):
            compound = compoundList[i]
            compound = removeStopWords(compound, stopWords)
            if getKey(compound) in compoundDictionary or getKey(compound) in unknownList or compound in passList:
                pass
            elif compound in [" ", "", "-", ":"]:
                pass
            else: # Name of compound not found in dictionary or lists
                print("Reading compound {} of {}".format(i+1, len(compoundList)))
                runAgain = True
                while(runAgain and getKey(compound) not in compoundDictionary):
                    runAgain = False
                    inputText = input("Enter the name of the following compound:\n{}\n$:".format(compound))
                    # PARSE INPUT
                    if inputText == INPUT_QUIT: # Quit
                        saveFiles()
                        sys.exit()
                    elif inputText == INPUT_QUIT_WITHOUT_SAVING: # Quit without saving
                        sys.exit()
                    elif inputText == INPUT_SAVE: # Save
                        saveFiles()
                        runAgain = True
                    elif inputText == INPUT_UNKNOWN: # Unknown compound
                        unknownList.append(getKey(compound))
                        history.append(i)
                    elif inputText == INPUT_STOP_WORDS:
                        ignored_words = compound.split(" ")
                        for word in ignored_words:
                            if word not in stopWords:
                                stopWords.append(word)
                                writeJson(stopWords, STOP_WORDS_FILE, indent=2)
                        history.append(i)
                    elif inputText == INPUT_UNDO:
                        if history == []:
                            print("Unable to undo")
                            runAgain = True
                        else:
                            oldIndex = history[-1]
                            oldNameKey = getKey(compoundList[oldIndex])
                            if oldNameKey in compoundDictionary:
                                oldValue = compoundDictionary[oldNameKey]
                                del compoundDictionary[oldNameKey]
                                print("Removed key from dictionary:\n{} : {}".format(oldNameKey, oldValue))
                            elif oldNameKey in unknownList:
                                del unknownList[index(oldNameKey)]
                                print("Removed {} from unknownList".format(oldNameKey))
                            saveFiles()
                            del history[-1]
                            i = oldIndex - 1
                    elif inputText[:len(INPUT_ADD_STOP_WORD)] == INPUT_ADD_STOP_WORD: # Add stop word
                        wordToAdd = inputText[len(INPUT_ADD_STOP_WORD)+1:]
                        if wordToAdd not in stopWords:
                            stopWords.append(wordToAdd)
                            writeJson(stopWords, STOP_WORDS_FILE, indent=2)
                            print("Added stop word {}".format(wordToAdd))
                            compound = removeStopWords(compound, stopWords)
                            runAgain = True
                    elif inputText == INPUT_PASS:
                        passList.append(compound)
                    elif inputText == "":
                        runAgain = True
                    else: # Normal input to add to dictionary
                        nameOfCompound = inputText
                        if inputText == INPUT_SAME:
                            nameOfCompound = compound
                        inputText = input("Add the following key to the dictionary? (Press n to cancel, ENTER to confirm):\n{} : {}\n$:".format(compound, nameOfCompound))
                        if inputText != "n":
                            compoundDictionary[getKey(compound)] = nameOfCompound
                            # Add the value to the dictionary with itself as the key
                            if autoAdd and getKey(nameOfCompound.lower()) not in compoundDictionary:
                                compoundDictionary[getKey(nameOfCompound.lower())] = nameOfCompound
                            history.append(i)
                            print("Added")
                        else:
                            runAgain = True
                if autoSave:
                    saveFiles()
        else: # END
            runAgain = True
            while(runAgain):
                runAgain = False
                inputText = input("Reached end of list. Save to dictionary? (Press n to cancel, ENTER to confirm, u to undo)\n$:")
                if inputText == "n":
                    inputText = input("Are you sure you want to quit without saving? (y/n)\n$:")
                    if inputText == "y":
                        sys.exit()
                    else:
                        saveFiles()
                elif inputText == INPUT_UNDO: # Undo
                    if history == []:
                        print("Unable to undo")
                        runAgain = True
                    else:
                        oldIndex = history[-1]
                        oldNameKey = getKey(compoundList[oldIndex])
                        if oldNameKey in compoundDictionary:
                            oldValue = compoundDictionary[oldNameKey]
                            del compoundDictionary[oldNameKey]
                            print("Removed key from dictionary:\n{} : {}".format(oldNameKey, oldValue))
                        elif oldNameKey in unknownList:
                            del unknownList[index(oldNameKey)]
                            print("Removed {} from unknownList".format(oldNameKey))
                        saveFiles()
                        del history[-1]
                        i = oldIndex - 1
                else: # Any other input
                    saveFiles()
                    sys.exit()
        i += 1

if __name__ == "__main__":
    structureList = loadStructures(STRUCTURES_FILE)
    # getCompressedDictionary(compoundDictionary, COMPRESSED_DICTIONARY_FILE)
    # parseAllDetails(structureList)
    # writeStructures(structureList, STRUCTURES_FILE)
    # compoundList = getCompoundList(structureList, useGetKey=False)
    # generateDictionary(compoundList)
    # printRecognizedCompounds(compoundList)
