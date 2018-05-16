#######################################
#***** X2K Web Genetic Algorithm *****#
#######################################

import os
import random
import PythonScripts.X2Kweb_API as xweb


# All possible X2K options
all_x2k_options = {
    'TF-target gene background database used for enrichment': [
        'ChEA 2015',
        'ENCODE 2015',
        'ChEA & ENCODE Consensus',
        'Transfac and Jaspar',
        'ChEA 2016',
        'ARCHS4 TFs Coexp',
        'CREEDS',
        'Enrichr Submissions TF-Gene Coocurrence',
    ],
    'kinase interactions to include': [#kea 2016
        'kea 2018',
        'ARCHS4',
        'iPTMnet',
        'NetworkIN',
        'Phospho.ELM',
        'Phosphopoint',
        'PhosphoPlus',
        'MINT',
    ],
    'enable_ppi': [
        'ppid',
        'Stelzl',
        'IntAct',
        'MINT',
        'BioGRID',
        'Biocarta',
        'BioPlex',
        'DIP',
        'huMAP',
        'InnateDB',
        'KEGG',
        'SNAVI',
        'iREF',
        'vidal',
        'BIND',
        'figeys',
        'HPRD',
    ],
    'max_number_of_interactions_per_article':  {"10":15, "01":50, "11":200, "00":1000000},
    'max_number_of_interactions_per_protein': {"10":50, "01":100, "11":200, "00":500},
    'min_network_size': {"10":1, "01":10, "11":50, "00":100},
    'min_number_of_articles_supporting_interaction': {"10":0, "01":1, "11":5, "00":10},
    'path_length': {"0":1, "1":2},
    'included organisms in the background database': {"10": "human", "01": "mouse", "11": "both", "00": "RESHUFFLE"},
}

def stringLength():
    string_length=0
    for key in all_x2k_options.keys():
        string_length += len(all_x2k_options[key])
    return  string_length
# stringLength()

def createPopulation(popSize, binaryStringLength):
    from random import choice
    populationinit = []
    for i in range(popSize):
        populationinit.append(''.join(choice(('0', '1')) for _ in range(binaryStringLength)) )
        print(populationinit[i])
    return populationinit

def reshuffle(x2k_options):
    new_options = x2k_options.copy()
    for param in new_options:
        selection = new_options[param]
        while selection == "RESHUFFLE":
            selection = random.choice( list(all_x2k_options[param].values()) )
        new_options[param] = selection
    return new_options

def translateDatabases(binaryString_segment, _dbs):
    selection = []
    for i, bit in enumerate(binaryString_segment):
        if bit == "1":
            selection.append(_dbs[i])
    return selection

def translateParameters(binaryString):
    x2k_options={}
    stringCount = 0
    # Database lists
    for db in ['TF-target gene background database used for enrichment','kinase interactions to include','enable_ppi']:
        dbList = all_x2k_options[db]
        selection = translateDatabases(binaryString[stringCount:stringCount + len(dbList)], dbList)
        x2k_options[db] = selection
        stringCount += len(selection)
    # Parameter dictionaries
    for param in ['max_number_of_interactions_per_article', \
                  'max_number_of_interactions_per_protein', \
                  'min_network_size',\
                    'min_number_of_articles_supporting_interaction',\
                   'path_length',\
                   'included organisms in the background database']:
        paramDict = all_x2k_options[param]
        bitLength = len(list(paramDict.keys())[0])
        bits = binaryString[stringCount:stringCount + bitLength]
        selection = paramDict[bits]

        x2k_options[param] = selection
        stringCount += bitLength
        # Reshuffle
        x2k_options = reshuffle(x2k_options)
    return x2k_options


"""
pop = createPopulation(10, stringLength() )
binaryString = pop[1]
"""


def calculateFitness(binaryString):
    kinase_file = 'Genetic_Algorithm/testgmt/' + os.listdir('Genetic_Algorithm/testgmt/')[0]
    # save_file = 'Genetic_Algorithm/output/x2k_out.csv'
    x2k_options = translateParameters(binaryString)

    x2k_output = xweb.run_X2K_allGenes(kinase_file, options=x2k_options, verbose=False, replaceNAs=False, outputValues='pvalue')

