import datetime
import os
import sys
import re
import json

gmt_id = {
    'CHEA.*': re.compile(
        r'^(?P<PMID>\d+)-(?P<Assay>.+?-.+?)-(?P<Cell_Line>.+?)_(?P<Species>.+)$'
    ),
    'CREEDS.*': re.compile(
        r'^(?P<Type>.+?)-(?P<GSID>GSE\d+)-sample-(?P<Sample>\d+)-CREEDS-.+?_(?P<Species>.+?)_(?P<UpOrDown>.+?)$'
    ),
    'Enrichr.*': re.compile(
        r'^Enrichr-Submissions_(?P<Species>.+?)$'
    ),
    'JASPAR-TRANSFAC.*': re.compile(
        r'^jaspar-transfac_(?P<Species>.+?)$'
    ),
}

def read_gmt_file(gmt_file):
    gmt_id_re = re.compile('.*')
    for k, v in gmt_id.items():
        if re.match(k, gmt_file):
            gmt_id_re = v
            break
    with open(gmt_file, 'r', encoding='utf-8') as fh:
        for line in fh:
            line_split = line.strip().split('\t')
            tf = line_split[0]
            tf_id = line_split[1]
            tf_id_parsed = {}
            try:
                tf_id_parsed = gmt_id_re.match(tf_id).groupdict()
            except Exception as e:
                print('WARN: Could not parse tf_id')
            genes = line_split[2:]
            yield {
                'tf': tf,
                'tf_id': tf_id,
                'genes': genes,
                'meta': tf_id_parsed,
            }

def get_timestamp():
    return datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')

def convert_gmt_row_to_background_row(gmt_row):
    for gene in gmt_row['genes']:
        yield [
            gmt_row['tf'],
            gmt_row['tf']+'-'+gmt_row['tf_id'],
            gene,
            'na',
            'na',
            'na',
            gmt_row.get('species', 'both'),
            get_timestamp(),
            json.dumps(gmt_row['meta']),
        ]

def write_background_file_with_gmt(background_file, gmt):
    with open(background_file, 'w', encoding='utf-8') as fh:
        ind = 0
        for gmt_row in gmt:
            for bg in convert_gmt_row_to_background_row(gmt_row):
                print(*([ind] + bg), sep='\t', end='\n', file=fh)
                ind += 1

for file in sys.argv[1:]:
    print('Processing', file)
    base, ext = os.path.splitext(file)
    write_background_file_with_gmt(
        base + '.tsv',
        read_gmt_file(file),
    )
