import pickle, json, os, sys, traceback
# import xml.etree.ElementTree as etree
from pathlib import Path
from pdb_crystal_database import loadStructures, writeStructures, getStructure, Structure
from misc_functions import loadJson, writeJson
from time import sleep

try:
    import requests
except ModuleNotFoundError:
    print("ERROR: Requests module not found. You must install the requests module in order to get structure info from the PDB."
    "For more information, see http://docs.python-requests.org/en/master/user/install/")
    sys.exit()

if not os.path.exists("Structures"):
    os.makedirs("Structures")

# Create Path objects for directories
INPUT_DIR = Path("Input/")
STRUCTURE_DIR = Path("Structures/")

WITHOUT_DETAILS_FILE = INPUT_DIR / "pdbs_without_details.json"
STRUCTURES_FILE = STRUCTURE_DIR / "structures.pkl" # The database file. Must be placed in proper location

structureList = []
pdbsWithoutDetails = []

# The number of pdb objects to download in one request. Kept at ~1000 to avoid going over the 8 KB request limit
# The script also saves structure info to disk every chunk
CHUNK_SIZE = 1000

# Will divide the list of all pdb ids into downloadable chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# Legacy XML code
# def loadPdb(pdbid): # ElementTree root
#     """Loads a pdbid as an xml file and returns the <root> branch"""
#     try:
#         response = requests.get("http://www.rcsb.org/pdb/rest/customReport.xml?pdbids="+pdbid+
#         "&customReportColumns=crystallizationMethod,crystallizationTempK,pdbxDetails,phValue,pmc,sequence,resolution&service=wsfile", timeout=10)
#     except (requests.Timeout, requests.exceptions.ConnectionError):
#         print("Request timeout on PDB: {}".format(pdbid))
#         return None
#     root = etree.fromstring(response.content)
#     if (root.find("record") != None):
#         return root
#     else:
#         print("\n-------------------- ERROR --------------------\nNo entry found with PDB ID "
#         + str(pdbid)+"\n-----------------------------------------------\n")
#         return None

def loadPdbs(pdbids):
    """ Loads pdbs as a dictionary object (from a JSON HTTP response) 
        pdbids: a list of PDB ID strings - ex. ["1STP", "2JEF", "1CDG"]  
    """

    # If you want to download more information, change this query
    # See https://data.rcsb.org/#gql-api
    # TODO Add sequences
    query = """
    {
        entries(entry_ids: [[PDB_LIST]]) {
            rcsb_id
            exptl_crystal_grow {
                method
                temp
                pdbx_details
                pH

            }
            rcsb_entry_info {
                resolution_combined
            }
            pubmed {
                rcsb_pubmed_central_id
            }
        }
    }
    """

    query = query.replace("[[PDB_LIST]]", str(pdbids)).replace("'", '"')
    # query = ''.join(query.split()) # Remove whitespace

    request = "https://data.rcsb.org/graphql?query=" + query

    # print("Request = " + request)

    response = requests.get(request)

    #TODO: error checking
    if (response.status_code != 200):
        print("Error with HTTP request")
        print(response.text)

    return json.loads(response.text)


def fetchStructures(pdbList, structureFile=STRUCTURES_FILE, onlyDetails=True, ignorePdbsWithoutDetails=True, ignoreCompletedPdbs=True): # list
    """Takes a list of pdbids and creates Structure objects for them, outputing them to structureFile
    If onlyDetails is True the function will only output Structures that have crystallization details
    If ignoreCompletedPdbs is True, then the function will read the Structure file and ignore Pdbs which already
        have Structures associated with them
    If ignoreCompletedPdbs is False, then the function will recheck and update every structure in the structure file
    If ignorePdbsWithoutDetails is True then the function will ignore Pdbs from WITHOUT_DETAILS_FILE
    saveFrequency is the number of pdbs downloaded before the files are saved
        increasing this number makes the script faster, but more unsaved data is lost when the script is exited
    """
    global structureList
    global pdbsWithoutDetails

    count = 0

    structureList = []
    pdbsWithoutDetails = []
    
    try:
        print("Loading PDBs without details...")
        pdbsWithoutDetails = loadJson(WITHOUT_DETAILS_FILE)
    except:
        print("File {} not found. One will be created.".format(WITHOUT_DETAILS_FILE))

    completedPdbList = [] # List of Pdbs already turned into structures or without details
    try:
        structureList = loadStructures(structureFile)
        if ignoreCompletedPdbs:
            for struc in structureList:
                completedPdbList.append(struc.pdbid)
    except FileNotFoundError:
        print("File {} not found. A new structure file will be created.".format(structureFile))

    ignoredPdbList = []

    if ignorePdbsWithoutDetails:
        ignoredPdbList.extend(pdbsWithoutDetails)
    # If ignoreCompletedPdbs=False, then completedPdbList will be empty
    ignoredPdbList.extend(completedPdbList)

    print("Removing completed pdbs (May take a minute)...")
    pdbList = list(set(pdbList)-set(ignoredPdbList))

    if pdbList == []:
        print("All PDBs have been ignored. No more PDBs need to be downloaded.")


    for chunk in chunks(pdbList, CHUNK_SIZE):
        sleep(0.10) # Add delay to avoid sending too many requests

        if count == 0:
            print("Downloading {} structure objects from the pdb".format(len(pdbList)))
        newCount = count + len(chunk)
        print("Downloading pdbs {} through {} of {}...".format(count, newCount, len(pdbList)))
        count = newCount

        print('Downloading data...')
        response = loadPdbs(chunk)
        print('Downloaded')

        #TODO - error checking
        if "data" not in response:
            print("ERROR in chunk")
            continue

        # For every PDB entry in the chunk
        for entry in response["data"]["entries"]:
            pdbid = entry["rcsb_id"]
            
            if (entry["exptl_crystal_grow"] == None):
                pdbsWithoutDetails.append(pdbid)
                continue

            crystalEntry = entry["exptl_crystal_grow"][0]

            details = crystalEntry["pdbx_details"]
            if details == None:
                pdbsWithoutDetails.append(pdbid)
                if onlyDetails:
                    continue
            
            try:
                pmcid = entry["pubmed"]["rcsb_pubmed_central_id"]
                pmcid = pmcid[3:] # Strip "PMC" from PMC1234567
            except:
                pmcid = None

            try:
                pH = float(crystalEntry["pH"])
            except:
                pH = None

            try:
                temperature = float(crystalEntry["temp"])
            except:
                temperature = None
            
            method = crystalEntry["method"]

            #TODO - sequences
            sequences = []

            try:
                resolution = float(entry["rcsb_entry_info"]["resolution_combined"][0])
            except:
                resultion = None

            structure = Structure(pdbid, pmcid, details, [], pH, temperature, method, sequences, resolution)
            # If the pdb already has a structure in the list, update it
            structureAlreadyInList = getStructure(structureList, structure.pdbid)
            if structureAlreadyInList != None:
                structureList.remove(structureAlreadyInList)
            structureList.append(structure)
            # print(str(structure))

        print('Saving data...')
        # Save structures
        writeStructures(structureList, structureFile)
        writeJson(pdbsWithoutDetails, WITHOUT_DETAILS_FILE)
        print('Saved')

    # Legacy XML
    # for pdbid in pdbList:
    #     sleep(0.05) # Add delay to avoid sending too many requests
    #     if count == 1:
    #         print("Downloading {} structure objects from the pdb".format(len(pdbList)))
    #     if count % 100 == 0:
    #         print("Loading pdb {} of {}...".format(count, len(pdbList)))
    #     root = loadPdb(pdbid)
    #     count += 1
    #     if (root != None):
    #         details = root.find("record/dimStructure.pdbxDetails").text
    #         if details == "null":
    #             details = None
    #             pdbsWithoutDetails.append(pdbid)
    #         if details != None or not onlyDetails:
    #             pmcid = root.find("record/dimStructure.pmc").text
    #             if (pmcid != "null"):
    #                 pmcid = pmcid[3:]
    #             else:
    #                 pmcid = None
    #             try:
    #                 pH = float(root.find("record/dimStructure.phValue").text)
    #             except ValueError:
    #                 pH = None
    #             try:
    #                 temperature = float(root.find("record/dimStructure.crystallizationTempK").text)
    #             except ValueError:
    #                 temperature = None
    #             method = root.find("record/dimStructure.crystallizationMethod").text
    #             if method == "null":
    #                 method = None
    #             sequences = []
    #             for tag in root.findall("record/dimEntity.sequence"):
    #                 sequences.append(tag.text)
    #             try:
    #                 resolution = float(root.find("record/dimStructure.resolution").text)
    #             except ValueError:
    #                 resolution = None
    #             structure = Structure(pdbid, pmcid, details, [], pH, temperature, method, sequences, resolution)
    #             # If the pdb already has a structure in the list, update it
    #             structureAlreadyInList = getStructure(structureList, structure.pdbid)
    #             if structureAlreadyInList != None:
    #                 structureList.remove(structureAlreadyInList)
    #             structureList.append(structure)
    #             if count % saveFrequency == 0:
    #                 # Save structures
    #                 writeStructures(structureList, structureFile)
    #                 writeJson(pdbsWithoutDetails, WITHOUT_DETAILS_FILE)

    writeStructures(structureList, structureFile)
    writeJson(pdbsWithoutDetails, WITHOUT_DETAILS_FILE)
    print("Done fetching Structures")
    return structureList

def getAllPdbs(filename=""): # list
    """Returns a list of every pdbid in the PDB
    # filename = an optional argument to also export the list as a json file"""
    print("Getting a list of all pdbs...")


    # Legacy XML
    # response = requests.get("https://www.rcsb.org/pdb/json/getCurrent").text
    # dictionary = json.loads(response)
    # allPdbs = dictionary["idList"]

    response = requests.get("https://data.rcsb.org/rest/v1/holdings/current/entry_ids")

    pdbids = json.loads(response.text)

    if filename != "":
        writeJson(pdbids, filename)

    return pdbids


if __name__ == "__main__":
    allPdbs = getAllPdbs()

    try:
        fetchStructures(allPdbs)
    except (Exception, KeyboardInterrupt) as e:
        # Print traceback
        print(traceback.format_exc())
        print("An error or keyboard interrupt was encountered. Saving downloaded data to disk...")
        writeStructures(structureList, STRUCTURES_FILE)
        writeJson(pdbsWithoutDetails, WITHOUT_DETAILS_FILE)
        print("Done")

    # result = loadPdbs(allPdbs[::1000])
    
    # writeJson(result, 'test.json', indent=4)
