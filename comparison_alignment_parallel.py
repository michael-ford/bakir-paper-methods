import pandas as pd
import dill as pickle
import sys
import os
import logging
from glob import glob
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from bakir.common import fasta_from_seq
from comparison_functions import (
    make_annotation_matches,
    group_annotations_by_gene,
    group_annotations_by_gene_immunanot,
    analyze_discordant_alleles,
    get_bakir_data,
    extract_immunannot_gene_features,
    load_bakir_data,
    load_skirt_data,
    compare_bakir_immunanot
)

logging.basicConfig(
    level=logging.WARNING,  # Set the logging level (DEBUG, INFO, WARNING, etc.)
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('comparison_alignment_parallel.log'),  # Save logs to the specified file
        logging.StreamHandler()  # Optional: also print to console
    ]
)

sys.path.insert(0, 'bakir/src/bakir/')

def load_data():
    with open('bakir/src/bakir/data/kir_db.pickle', 'rb') as f:
        anno = pickle.load(f)
    return anno[0]

def prepare_dfs(assembly_files: list):
    dfs = {}
    dfs_s = {}
    
    for assembly in assembly_files:
        bakir = load_bakir_data(assembly)
        prefix = os.path.basename(assembly).split('.f1')[0]

        immunanot = extract_immunannot_gene_features(f'HPRC-Immunanot-annotations/{prefix}/{prefix}.gtf.gz')
        dfs[prefix] = compare_bakir_immunanot(bakir, immunanot)

        skirt_file = glob(f'HPRC-Skirt-annotations//{prefix}/{prefix}*allele.csv')[0]
        skirt = load_skirt_data(skirt_file)
        dfs_s[prefix] = compare_bakir_immunanot(bakir, skirt)
    
    return dfs, dfs_s

def process_sample(k, dfs, dfs_s, db):
    logging.info(f"---> Sample {k}")
    all_data = []

    diff_cnv_i = dfs[k][dfs[k]['kir_annots copies'] != dfs[k]['imm_annots copies']]
    diff_a_i = dfs[k][dfs[k]['num allele call diffs'] > 0]

    diff_cnv_s = dfs_s[k][dfs_s[k]['kir_annots copies'] != dfs_s[k]['imm_annots copies']]
    diff_a_s = dfs_s[k][dfs_s[k]['num allele call diffs'] > 0]

    genes = set(list(diff_cnv_i['gene']) + list(diff_a_i['gene']) + list(diff_cnv_s['gene']) + list(diff_a_s['gene']))
    
    if genes:
        for g in genes:
            logging.info(f"----> Gene: {g}")
            sample = k.split('.')[0]
            haplo = k
            ka_geno = get_bakir_data(sample, haplo)
            ka_geno = group_annotations_by_gene(ka_geno)[g]

            alt_i = group_annotations_by_gene_immunanot(
                extract_immunannot_gene_features(f'HPRC-Immunanot-annotations/{haplo}/{haplo}.gtf.gz')
            )[g]
            alt_s = group_annotations_by_gene_immunanot(
                load_skirt_data(glob(f'HPRC-Skirt-annotations/{haplo}/{haplo}*allele.csv')[0])
            )[g]

            match_indexes, anno_matches_i, _, _ = make_annotation_matches(ka_geno, alt_i)
            match_indexes, anno_matches_s, _, _ = make_annotation_matches(ka_geno, alt_s)

            for ka_anno, im_anno in anno_matches_i:
                try:
                    ka_allele = ka_anno['closest allele'].split('*')[1][:3]
                    im_allele = im_anno['template_allele'].split('*')[1][:3]

                    if anno_matches_s:
                        s_anno = [x[1] for x in anno_matches_s if ka_anno == x[0]][0]
                        s_allele = s_anno['template_allele'].split('*')[1][:3]
                    else:
                        s_anno = []
                        s_allele = None

                    if (ka_allele == im_allele and im_allele == s_allele) and s_allele:
                        continue

                    data = OrderedDict([
                        ('sample', sample),
                        ('haplotype', k.split('.')[1]),
                        ('contig', ka_anno['reference']),
                        ('BAKIR allele', ka_anno['closest allele']),
                        ('BAKIR positions', (ka_anno['start'], ka_anno['end']))
                    ])

                    closest_alleles_i, _, _ = analyze_discordant_alleles(ka_anno, im_anno, db=db, display_table=False)
                    closest_alleles_i = closest_alleles_i.set_index('closest allele')

                    closest_alleles_s_f = False
                    if '-' in s_anno['gene_name']:
                        closest_alleles_s = dict([(x, None) for x in ['common mut', 'new mut', 'missing mut', 'missing functional mut',
                                                                       'new functional mut', 'jaccard functional distance', 'jaccard distance',
                                                                       'MSA edit distance']])
                        closest_alleles_s['closest allele'] = s_anno['template_allele']
                        closest_alleles_s = pd.DataFrame([closest_alleles_s])
                        closest_alleles_s = closest_alleles_s.set_index('closest allele')
                        closest_alleles_s_f = True
                    elif s_allele and (s_allele != im_allele):
                        closest_alleles_s, _, _ = analyze_discordant_alleles(ka_anno, s_anno, db=db, display_table=False)
                        closest_alleles_s = closest_alleles_s.set_index('closest allele')
                        closest_alleles_s_f = True

                    if closest_alleles_s_f:
                        closest_alleles = pd.concat([closest_alleles_i, closest_alleles_s])
                        closest_alleles = closest_alleles.reset_index().drop_duplicates(subset='closest allele').set_index('closest allele')
                    else:
                        closest_alleles = closest_alleles_i

                    data[f'BAKIR global edit distance'] = closest_alleles.loc[ka_anno['closest allele']]['MSA edit distance']
                    data[f'BAKIR core jaccard distance'] = closest_alleles.loc[ka_anno['closest allele']]['jaccard functional distance']
                    data[f'BAKIR missing core variants'] = closest_alleles.loc[ka_anno['closest allele']]['missing functional mut']
                    data[f'BAKIR new core variants'] = closest_alleles.loc[ka_anno['closest allele']]['new functional mut']

                    for l, an in zip(['Immunannot', 'SKIRT'], [im_anno, s_anno]):
                        data[f'{l} allele'] = an['template_allele'] if an else None
                        data[f'{l} positions'] = (an['start'], an['end']) if an else None
                        data[f'{l} global edit distance'] = closest_alleles.loc[an['template_allele']]['MSA edit distance'] if an else None
                        data[f'{l} core jaccard distance'] = closest_alleles.loc[an['template_allele']]['jaccard functional distance'] if an else None
                        data[f'{l} missing core variants'] = closest_alleles.loc[an['template_allele']]['missing functional mut'] if an else None
                        data[f'{l} new core variants'] = closest_alleles.loc[an['template_allele']]['new functional mut'] if an else None

                    all_data.append(data)

                except Exception as e:
                    logging.exception(f"Exception in {k}, gene {g}", exc_info=e)
    
    return all_data

def main(max_workers: int = 10):
    db = load_data()
    assembly_files = glob("HPRC-assemblies-annotations/*/*.*.f1_assembly_v2.yaml")
    dfs, dfs_s = prepare_dfs(assembly_files)

    all_data = []
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_sample, k, dfs, dfs_s, db): k for k in dfs}
        for future in as_completed(futures):
            k = futures[future]
            try:
                data = future.result()
                all_data.extend(data)
            except Exception as e:
                logging.exception(f"Exception during processing sample {k}", exc_info=e)

    pd.DataFrame(all_data).to_csv('all_discordant_comparison.csv', index=False)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run parallel analysis of discordant alleles.")
    parser.add_argument('--max-workers', type=int, default=10, help="Maximum number of parallel processes.")
    
    args = parser.parse_args()
    
    main(max_workers=args.max_workers)
