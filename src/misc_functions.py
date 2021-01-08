import json

def loadJson(jsonFile, object_pairs_hook=None):
    """Loads a json file, returns whatever python object is in the file
    Optional arguments are the same as in json.load"""
    with open(jsonFile, "r") as inFile:
        return json.load(inFile, object_pairs_hook=object_pairs_hook)

def writeJson(object, filename, indent=None, sort_keys=False):
    """Writes a python object to a json file
    Optional arguments are the same as in json.dumps"""
    with open(filename, "w") as outFile:
        outFile.write(json.dumps(object, indent=indent, sort_keys=sort_keys))

def printList(list, delimiter=", "): # string
    """Takes a list and returns a string deliminated by a specific substring"""
    string = ""
    count = 1
    for element in list:
        string += element
        if (count != len(list)):
            string += delimiter
        count += 1
    return string

def fileToList(filename):
    """Returns a list of strings, with each element being a line from the file"""
    list = []
    with open(filename, "r") as fileOfList:
        for element in fileOfList:
            list.append(element.rstrip())
    return list

def listToFile(list, filename): # void
    """Writes a list to a text file, with each element being a line in the file"""
    with open(filename, "w") as outFile:
        for element in list:
            outFile.write(str(element)+"\n")

def getKey(compound): # string
    """Takes a string name of a compound and removes certain characters to make a simplified dictionary key"""
    charactersToRemove = [" ", "(", ")", "-", ",", "/", "*", "_"]
    for char in charactersToRemove:
        compound = compound.replace(char, "")
    return compound.lower()
