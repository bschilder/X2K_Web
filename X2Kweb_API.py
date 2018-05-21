

"""
1. Generate random population 
2. Fix and mutate 
3. 


"""


import pandas as pd
import traceback
import random
import numpy as np
import http.client
import json
import os
import sys
import time



def timeit(func, *kargs, **kwargs):
    start = time.clock()
    func(*kargs, **kwargs)
    return time.clock() - start

class Pick:
    def __init__(self, n_picks, choices):
        self.n_picks = n_picks
        self.choices = choices

    def next(self):
        return random.sample(self.choices, k=self.n_picks)

# Randomly generate three options
def test_Pick():
    L = list(range(5))
    for n in range(5):
        l = []
        for val in Pick(3, L).next():
            assert val in L
            assert val not in l, "val (%d) was in l (%s)" % (val, l)
            l.append(val)

def enrichr_to_input_genes_target_kinase(fn):
    with open(fn, 'r', encoding='utf8') as fh:
        for line in fh:
            line_split = line.split('\t')
            kinase = line_split[0].split('_')[0]
            genes = list(filter(None, [
                gene_val.split(',')[0].strip()
                for gene_val in line_split[2:]
            ]))
            yield (kinase, genes)



def test_enrichr_to_input_genes_target_kinase():
    global test_kinase, test_input_genes
    for kinase, genes in enrichr_to_input_genes_target_kinase(test_kinase_file):
        test_kinase = kinase
        print(kinase)
        test_input_genes = genes
        print(genes)
        break


def produce_random_option(options_grid):
    return {
        option_name: option_value.next() if isinstance(option_value, Pick) else option_value
        for option_name, option_value in options_grid.items()
    }

def test_produce_random_option():
    global test_options
    test_options = produce_random_option(options_grid)
    print(test_options)


def prepare_options_for_x2k(input_genes, options):
    options = options.copy()

    # Add input_genes
    options['text-genes'] = input_genes

    # Convert ppi into enable flags
    for ppi in options['enable_ppi']:
        options['enable_' + ppi] = 'true'
    del options['enable_ppi']

    # Convert any lists
    return {
        k: '\n'.join(v) if type(v) == list else str(v)
        for k, v in options.items()
    }

def test_prepare_options_for_x2k():
    global test_x2k_options
    test_x2k_options = prepare_options_for_x2k(
        test_input_genes, test_options
    )
    print(
        test_x2k_options
    )


def create_x2k_payload(x2k_options):
    boundary = 'WebKitFormBoundary3ruqqYyh1YAF9JXu'
    payload = '------' + boundary + '\r\n' \
        + ('------' + boundary + '\r\n').join(
            'Content-Disposition: form-data; name="' \
            + key + '"\r\n\r\n' \
            + value + '\r\n'
            for key, value in x2k_options.items()
        ) + '------' + boundary + '--'
    headers = {
        'content-type': 'multipart/form-data; boundary=----' + boundary,
        'cache-control': 'no-cache',
    }
    return (headers, payload)

def test_create_x2k_payload():
    print(*create_x2k_payload(test_x2k_options), sep='\n')


def query_x2k(x2k_options):
    conn = http.client.HTTPConnection("localhost:8080", timeout=20)
    # conn = http.client.HTTPConnection("amp.pharm.mssm.edu", timeout=15)

    headers, payload = create_x2k_payload(x2k_options)
    conn.request("POST", "/X2K/api", payload, headers)

    res = conn.getresponse()
    if res.status != 200:
        raise Exception("Server response was %d" % (res.status))

    data = res.read().decode('utf-8')

    conn.close()

    return (res, data)

def test_query_x2k():
    global test_x2k_response
    res, test_x2k_response = query_x2k(
        test_x2k_options
    )
    print(
        res,
        test_x2k_response,
        sep='\n',
    )


def parse_x2k_response_kinases(x2k_results, significant=0.05):
    x2k_results_parsed = {
        key: json.loads(value) if key != 'input' else value
        for key, value in json.loads(x2k_results).items()
    }
    kinases = {
        kinase['name']: kinase
        for kinase in x2k_results_parsed['KEA']['kinases']
    }
    n_g2n_nodes = len([
        tf
        for tf in x2k_results_parsed['G2N']['network']['nodes']
        if tf['pvalue'] < significant
    ])
    n_sig_kinases = len([
        kinase for kinase in kinases.values()
        if kinase['pvalue'] < significant
    ])
    return {
        'kinases': kinases,
        'n_g2n_nodes': n_g2n_nodes,
        'n_sig_kinases': n_sig_kinases,
    }

def test_parse_x2k_response_kinases():
    global test_x2k_results
    test_x2k_results = parse_x2k_response_kinases(test_x2k_response)
    print(
        test_x2k_results
    )

def evaluate_genelist_on_x2k(input_genes, options):
    res, data = query_x2k(
        prepare_options_for_x2k(
            input_genes, options
        )
    )
    return parse_x2k_response_kinases(data)

def test_evaluate_genelist_on_x2k():
    print(
        evaluate_genelist_on_x2k(
            test_input_genes,
            test_options,
        ),
        sep='\n',
    )





"""
results = x2k_results
kinase='WEE1'
"""

"""
def replace_NAs_with_random_rank(DF):
    for col in DF.iloc[:,1:]:
        print(col)
        maxRank = int(max(DF[col].dropna()))
        remainingRanks = list(range(maxRank+1,len(DF)+1))
        shuffle(remainingRanks)
        DF[col] = DF[col].apply(choose_random_rank).astype(int)
    #DF.to_csv(outputName, sep='\t', header=True, index=None, na_rep='NA')
    return DF 

def choose_random_rank(input):
    if math.isnan(input):
        newRank = np.random.choice(remainingRanks, 1, replace=False)[0]
    else:
        newRank = input
    return newRank
"""

def return_all_kinase_ranks(results, missing_score=1):
    sortedResults = pd.DataFrame(results['kinases']).transpose()
    sortedResults.head()
    sortedResults = sortedResults.sort_values(by='pvalue', ascending=True).reset_index()
    sortedResults['Rank'] = range(0,len(sortedResults))
    return sortedResults[['name','Rank']]


# Calculate z-scores from each
def return_all_kinase_pvalues(results, missing_score=1):
    sortedResults = pd.DataFrame(results['kinases']).transpose()
    sortedResults.head()
    sortedResults = sortedResults.sort_values(by='pvalue', ascending=True).reset_index()
    #sortedResults['-log(pvalue)'] = -np.log(pd.to_numeric(sortedResults['pvalue']))
    pvalue_DF = sortedResults[['name', 'pvalue']]
    pvalue_DF.columns = ['Kinase','pvalue']
    return pvalue_DF

def test_return_all_kinase_ranks():
    assert return_all_kinase_ranks({}, -1) == "hi"


def score_kinase(kinase, results, missing_score=1):
    try:
        val = -np.log(results['kinases'][kinase]['pvalue'])
        #zscore = results['kinases'][kinase]['zscore'] # zscore across kinases, not across experiments?
        if val < 1:
            raise Exception("How? " + results['kinases'][kinase]['pvalue'])
        return val
    except:
        return missing_score

def test_score_kinase():
    global test_kinase_score
    kinase = list(test_x2k_results['kinases'].keys())[0]
    test_kinase_score = score_kinase(
        kinase,
        test_x2k_results,
    )
    print(
        kinase,
        test_kinase_score,
    )


def options_to_binary_array(options, options_grid):
    keys = sorted(options_grid.keys())
    key_choices = {option: definition.choices
                   for option, definition in options_grid.items()
                   if isinstance(definition, Pick)}
    choices_keys = [
        '%s.%s' % (option, choice)
        for option, choices in key_choices.items()
        for choice in choices
    ]
    output = []
    for key, choices in zip(keys, map(key_choices.get, keys)):
        if choices is not None:
            output += np.isin(choices, options[key]).astype(np.int8).tolist()
    return pd.DataFrame(output, index=choices_keys).T

def test_options_to_binary_array():
    global test_options_binary_array
    test_options_binary_array = options_to_binary_array(
        test_options,
        options_grid,
    )
    print(
        test_options,
        options_grid,
        test_options_binary_array.T,
        sep='\n',
    )


def calculate_summary_results(results):
    values = results \
        .drop(['score'], axis=1) \
        .dropna()
    unweighted = values \
        .sum(axis=0)
    weighted = values \
        .multiply(results['score'], axis=0) \
        .sum(axis=0)
    return (weighted / unweighted) \
        .to_frame() \
        .T

def test_calculate_sumnmary_results():
    global test_summary_result
    test_summary_result = test_options_binary_array.copy()
    test_summary_result['score'] = test_kinase_score
    print(test_summary_result)
    print(
        test_summary_result,
        calculate_summary_results(test_summary_result),
        sep='\n',
    )


def first(arr):
    for a in arr:
        return a

def convert_result_to_list(result, header=False):
    if header:
        return first(result.iterrows())[1].index.tolist()
    else:
        return first(result.iterrows())[1].tolist()

def test_convert_result_to_list():
    print(
        convert_result_to_list(test_summary_result, True),
        convert_result_to_list(test_summary_result, False),
        sep='\n',
    )

""" TESTS
print(timeit(test_Pick))
test_kinase_file = 'GEO.txt'
print(timeit(test_enrichr_to_input_genes_target_kinase))
print(timeit(test_produce_random_option))
print(timeit(test_prepare_options_for_x2k))
print(timeit(test_create_x2k_payload))
print(timeit(test_query_x2k))
print(timeit(test_parse_x2k_response_kinases))
print(timeit(test_evaluate_genelist_on_x2k))
print(timeit(test_score_kinase))
print(timeit(test_options_to_binary_array))
print(timeit(test_calculate_sumnmary_results))
print(timeit(test_convert_result_to_list))
"""

# Use best options
best_options = {
    'TF-target gene background database used for enrichment': [
        # 'ChEA 2015',
        # 'ENCODE 2015',
        'ChEA & ENCODE Consensus',
        # 'Transfac and Jaspar',
        # 'ChEA 2016',
        # 'ARCHS4 TFs Coexp',
        # 'CREEDS',
        # 'Enrichr Submissions TF-Gene Coocurrence',
    ],
    'kinase interactions to include': [#kea 2016
        'kea 2018',
        # 'ARCHS4',
        # 'iPTMnet',
        # 'NetworkIN',
        # 'Phospho.ELM',
        # 'Phosphopoint',
        # 'PhosphoPlus',
        # 'MINT',
    ],
    'enable_ppi': [
        'ppid',
        'Stelzl',
        'IntAct',
        'MINT',
        'BioGRID',
        # 'Biocarta',
        # 'BioPlex',
        # 'DIP',
        # 'huMAP',
        # 'InnateDB',
        # 'KEGG',
        # 'SNAVI',
        # 'iREF',
        # 'vidal',
        # 'BIND',
        # 'figeys',
        # 'HPRD',
    ],
    'max_number_of_interactions_per_article': 1000000,
    'max_number_of_interactions_per_protein': 200,
    'min_network_size': 10,
    'min_number_of_articles_supporting_interaction': 0,
    'path_length': 2,
    'included organisms in the background database': 'both',
}

kinase_file = 'GEO.txt'
save_file = 'kinaseRank_X2K.txt'

# Only use GEO-contained kinases
# kinases = set(list(zip(*list(enrichr_to_input_genes_target_kinase('GEO.txt'))))[0])
# kinase_genes = [
#     (kinase, genes)
#     for kinase, genes in list(enrichr_to_input_genes_target_kinase(kinase_file))
#     if kinase in kinases
# ]

# For both X2K and KEA
def absent_perturbed_kinases(DF):
    perturbedKinases = set([x.split('_')[0] for x in DF.columns[1:]])
    absentPerturbed=[]
    for pert in perturbedKinases:
        syns = synDict[pert]
        overlap = set(syns).intersection(set(DF['Kinase']))
        if len(overlap)==0:
            absentPerturbed.append(pert)
    return absentPerturbed


"""
def run_X2K(kinase_file, save_file, options=best_options, verbose=False):
    save_file_exists = os.path.isfile(save_file)
    df = pd.read_csv(save_file, index_col=False).drop(['n_g2n_nodes', 'n_sig_kinase', 'kinase'],
                                                      axis=1) if save_file_exists else pd.DataFrame()
    kinase_genes = list(enrichr_to_input_genes_target_kinase(kinase_file))

    for kinase, genes in kinase_genes:
        try:
            last = time.time()
            # options = best_options
            n_g2n_nodes = float('nan')
            n_sig_kinase = float('nan')
            score = float('nan')
            try:
                x2k_results = evaluate_genelist_on_x2k(genes, options)
                score = score_kinase(kinase, x2k_results)
                n_g2n_nodes = x2k_results['n_g2n_nodes']
                n_sig_kinase = x2k_results['n_sig_kinases']
            except ConnectionRefusedError:
                time.sleep(20)
                continue
            except KeyboardInterrupt:
                break
            except Exception:
                traceback.print_exc(file=sys.stderr)
                score = float('nan')

            # Save result to file
            with open(save_file, 'a') as fh:
                if not save_file_exists:
                    print('score', 'n_g2n_nodes', 'n_sig_kinase', 'kinase',
                        sep=',', end='\n',
                        file=fh,
                    )
                    save_file_exists = True
                print(
                    score, n_g2n_nodes, n_sig_kinase, kinase,
                    sep=',', end='\n',
                    file=fh,
                )

            if verbose == True:
                # Calculate and print running summary
                print(kinase)
                print(
                    'Last Kinase: %s' % (kinase),
                    'Last score: %s' % (float(score)),
                    'Time since last result: %0.3f' % (time.time() - last),
                    sep='\n',
                )
                last = time.time()
        except KeyboardInterrupt:
            break
"""

#cd Kinase_Enrichment_Comparisons

def run_X2K_allGenes(kinase_file, save_file='No', options=best_options, verbose=False, outputValues='pvalue'):
    finalDF = pd.DataFrame(columns=['Kinase'])
    try:
        os.remove(save_file)
    except OSError:
        pass
    kinase_genes = list(enrichr_to_input_genes_target_kinase(kinase_file))
    experiment_names = pd.read_table(kinase_file, header=None).iloc[:,0].tolist()
    index = 0
    for kinase, genes in kinase_genes:
        print(kinase)
        try:
            last = time.time()
            score = float('nan')
            try:
                x2k_results = evaluate_genelist_on_x2k(genes, options)
                if outputValues=='pvalue':
                    kinase_values = return_all_kinase_pvalues(x2k_results, missing_score=1)
                elif outputValues=='pvalue ranks':
                    kinase_values = return_all_kinase_ranks(x2k_results, missing_score=1)
                    kinase_values.columns = ['Kinase', experiment_names[index]]
                finalDF = finalDF.merge(kinase_values, on='Kinase', how='outer')

            except ConnectionRefusedError:
                time.sleep(20)
                continue
            except KeyboardInterrupt:
                break
            except Exception:
                traceback.print_exc(file=sys.stderr)
                score = float('nan')

            if verbose == True:
                print(experiment_names[index])
                # Calculate and print running summary
                print(kinase)
                # print('Last Kinase: %s' % (kinase),
                #     'Last score: %s' % (float(score)),
                #     'Time since last result: %0.3f' % (time.time() - last),
                #     sep='\n')

        except KeyboardInterrupt:
            break
        index+=1
    # Save result to file
    if save_file!='No':
        finalDF.to_csv(save_file, sep='\t', header=True, index=None, na_rep='NA')
    return finalDF



