import pandas as pd
from typing import List, Dict, Tuple, Set, Any
import yaml
from collections import defaultdict
import difflib
import gzip
import logging
import numpy as np

# Define a format to include the file name and line number
log_format = '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'

# Configure your logger
logging.basicConfig(level=logging.INFO, format=log_format)
logger = logging.getLogger("comparison_functions")



def load_skirt_data(file_path: str) -> List[Dict[str, str]]:
    """ Load the SKIRT data from a CSV file and return it as a list of dictionaries. 
    Args:
        file_path (str): The path to the CSV file.
    Returns:
        List[Dict[str, str]]: The SKIRT data as a list of dictionaries.
    """
    skirt = pd.read_csv(file_path)
    skirt = skirt.rename({'kir_allele': 'template_allele', 'target_start': 'start', 'target_end': 'end', 'target_name': 'contig'}, axis='columns')
    skirt['gene_name'] = skirt['template_allele'].apply(lambda x: x.split('*')[0])
    skirt = skirt.to_dict('records')
    return skirt

def extract_immunannot_gene_features(file_path: str) -> List[Dict[str, str]]:
    """
    Extracts features of all 'gene' lines from a GTF-like file compressed with gzip (.gtf.gz).

    Args:
        file_path (str): Path to the gzipped file containing the data.

    Returns:
        List[Dict[str, str]]: A list of dictionaries, where each dictionary contains attributes of a gene.
    """
    gene_features = []

    with gzip.open(file_path, 'rt') as file:  # Note the 'rt' mode for text mode reading
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue  # Skip comments and empty lines

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # Ensure the line has enough parts to process

            feature_type = parts[2]
            if feature_type == 'gene':
                attributes = parts[8]
                attributes_dict = {attr.strip().split(' ')[0]: attr.strip().split(' ')[1].strip('"') for attr in attributes.split(';') if attr}
                gene_features.append({
                    'contig': parts[0],
                    'source': parts[1],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6]
                })
                gene_features[-1].update(attributes_dict)

    return gene_features


def load_bakir_data(yaml_path: str) -> List[Any]:
    """
    Loads a YAML file, modifies each 'start' value by +1, and returns the modified data.
    
    Args:
        yaml_path (str): Path to the YAML file.
    
    Returns:
        List[Any]: The modified YAML data.
    """
    # Load the YAML data using UnsafeLoader
    with open(yaml_path, 'r') as file:
        data = yaml.load(file, Loader=yaml.UnsafeLoader)
    
    # Iterate through the list of ordered dictionaries and modify 'start' value
    for item in data:
        if 'start' in item:
            item['start'] += 1
    
    return data


def compare_bakir_immunanot(bakir_data: List[Dict], alt_annotation_data: List[Dict], alt_annots_label: str = "imm_annots") -> pd.DataFrame:
    """
    Compare KIR gene annotations from bakir and immunanot.
    """
    # Initialize an empty DataFrame for all discrepancies
    all_discrepancies_df = pd.DataFrame(columns=['gene', 'kir_annots copies', f'{alt_annots_label} copies', 'num position diffs', 'num strand diffs'])

    # Group annotations by gene for both tools
    kir_annotations = group_annotations_by_gene(bakir_data)
    imm_annotations = group_annotations_by_gene_immunanot(alt_annotation_data)

    for gene, kir_annots in kir_annotations.items():
        alt_annots = imm_annotations.get(gene, [])
        gene_discrepancies_df = compare_gene_annotations(gene, kir_annots, alt_annots, alt_annots_label)
        all_discrepancies_df = pd.concat([all_discrepancies_df, gene_discrepancies_df], ignore_index=True)

    return all_discrepancies_df

def group_annotations_by_gene(annotations: List[Dict]) -> Dict[str, List[Dict]]:
    """
    Group bakir annotations by gene.
    
    Args:
        annotations (List[Dict]): Annotations to group.
    
    Returns:
        Dict[str, List[Dict]]: Grouped annotations by gene.
    """
    grouped = defaultdict(list)
    for annot in annotations:
        grouped[annot['gene']].append(annot)
    return grouped

def group_annotations_by_gene_immunanot(annotations: List[Dict]) -> Dict[str, List[Dict]]:
    """
    Group immunanot annotations by gene.
    
    Args:
        annotations (List[Dict]): Annotations to group.
    
    Returns:
        Dict[str, List[Dict]]: Grouped annotations by gene.
    """
    grouped = defaultdict(list)
    for annot in annotations:
        gene_name = annot['gene_name']
        grouped[gene_name].append(annot)
    return grouped

def report_discrepancies(discrepancies: List[str]) -> None:
    """
    Report discrepancies found during comparison.
    
    Args:
        discrepancies (List[str]): Descriptions of discrepancies.
    """
    if discrepancies:
        print("Discrepancies found:")
        for discrepancy in discrepancies:
            print(discrepancy)
    else:
        print("No discrepancies found.")


def calculate_match_scores(kir_annots: List[Dict], alt_annots: List[Dict]) -> Dict[Tuple[int, int], int]:
    """
    Calculate match scores between bakir and immunanot annotations for a gene, returning a dictionary.
    
    Args:
        kir_annots (List[Dict]): Annotations from bakir.
        alt_annots (List[Dict]): Annotations from immunanot.
    
    Returns:
        Dict[Tuple[int, int], int]: A dictionary where keys are tuples of indexes (bakir, immunanot)
                                    and values are the match scores.
    """
    match_scores = {}
    for i, kir_annot in enumerate(kir_annots):
        for j, imm_annot in enumerate(alt_annots):
            start_diff = abs(int(kir_annot['start']) - int(imm_annot['start']))
            end_diff = abs(int(kir_annot['end']) - int(imm_annot['end']))
            score = start_diff + end_diff
            match_scores[(i, j)] = score
    return match_scores

def find_best_matches(match_scores: Dict[Tuple[int, int], int], kir_len: int, imm_len: int) -> Tuple[List[Tuple[int, int]], List[int], List[int]]:
    """
    Find the best matches and identify extra copies based on match scores, where match_scores is now a dictionary.
    
    Args:
        match_scores (Dict[Tuple[int, int], int]): Match scores between annotations, with keys as (bakir index, immunanot index) and values as scores.
        kir_len (int): Number of bakir annotations.
        imm_len (int): Number of immunanot annotations.
    
    Returns:
        Tuple[List[Tuple[int, int]], List[int], List[int]]: A tuple containing a list of best matches (pairs of indexes),
                                                             and lists of indexes for extra bakir and immunanot annotations.
    """
    # Convert match_scores dictionary to a sorted list of tuples by score
    sorted_scores = sorted(match_scores.items(), key=lambda x: x[1])
    
    best_matches, kir_matched, imm_matched = [], set(), set()
    for (i, j), _ in sorted_scores:
        if i not in kir_matched and j not in imm_matched:
            best_matches.append((i, j))
            kir_matched.add(i)
            imm_matched.add(j)
    
    # Identify extra copies
    extra_kir = [i for i in range(kir_len) if i not in kir_matched]
    extra_imm = [j for j in range(imm_len) if j not in imm_matched]
    
    return best_matches, extra_kir, extra_imm


def make_annotation_matches(kir_annots: List[Dict], alt_annots: List[Dict]):
    match_scores = calculate_match_scores(kir_annots, alt_annots)
    best_matches, extra_kir, extra_imm = find_best_matches(match_scores, len(kir_annots), len(alt_annots))

    return best_matches, [(kir_annots[i], alt_annots[j]) for i, j in best_matches], extra_kir, extra_imm

def compare_gene_annotations(gene: str, kir_annots: List[Dict], alt_annots: List[Dict], alt_annots_label: str = "imm_annots") -> pd.DataFrame:
    """
    Compare annotations for a specific gene between bakir and alt_annots_label based on closest start and end positions.
    """
    match_scores = calculate_match_scores(kir_annots, alt_annots)
    best_matches, extra_kir, extra_imm = find_best_matches(match_scores, len(kir_annots), len(alt_annots))
    # Initialize a list to collect row data
    discrepancies_data = []


    discrepancies_data.append({
        'gene': gene,
        'kir_annots copies': len(kir_annots),
        f'{alt_annots_label} copies': len(alt_annots)})

    position_diffs = 0
    position_diffs_values = []
    strand_diff = 0
    allele_diff = 0
    allele_diff_values = []
    for i, j in best_matches:
        score = match_scores[(i, j)]  # Access the score using the dictionary
        start_diff = kir_annots[i]['start'] - int(alt_annots[j]['start'])
        end_diff = kir_annots[i]['end'] - int(alt_annots[j]['end'])
        if start_diff != 0:
            logger.debug(f"Start position mismatch for {gene}: bakir copy {i},  alt_annots_label copy {j}: bakir start = {kir_annots[i]['start']}, alt_annots_label start = {alt_annots[j]['start']}.")
        if end_diff != 0:
            logger.debug(f"End position mismatch for {gene}: bakir copy {i},  alt_annots_label copy {j}: bakir end = {kir_annots[i]['end']}, alt_annots_label end = {alt_annots[j]['end']}.")
        
        if kir_annots[i]['is reverse'] != (alt_annots[j]['strand'] == '-'):
            logger.debug(f"Strand mismatch for {gene}: bakir copy {i},  alt_annots_label copy {j}: bakir strand = {'-'  if kir_annots[i]['is reverse'] else '+'}, alt_annots_label strand = {alt_annots[j]['strand']}.")
            strand_diff += 1
        
        if kir_annots[i]['closest allele'].split('*')[1][:3] != alt_annots[j]['template_allele'].split('*')[1][:3]:
            logger.debug(f"Allele mismatch for {gene}: bakir copy {i},  alt_annots_label copy {j}: bakir allele = {kir_annots[i]['closest allele']}, alt_annots_label allele = {alt_annots[j]['template_allele']}.")
            allele_diff += 1
            allele_diff_values.append((kir_annots[i]['closest allele'].split('*')[1], alt_annots[j]['template_allele'].split('*')[1]))


        # Add to df reporting discrepancies
        is_position_diff = (start_diff != 0 or end_diff != 0)
        position_diffs += int(is_position_diff)
        if is_position_diff:
            position_diffs_values.append((start_diff, end_diff))

    discrepancies_data[0].update({
        'num position diffs': position_diffs,
        'position diffs': position_diffs_values,
        'num strand diffs': int(strand_diff),  # 1 if there's a strand mismatch, else 0
        'num allele call diffs': int(allele_diff),
        'allele diffs': allele_diff_values
    })


    # Create DataFrame from collected row data
    discrepancies_df = pd.DataFrame(discrepancies_data)

    return discrepancies_df



import typing as t

def parse_bed_line(line: str) -> t.Dict[str, t.Any]:
    """Parse a single BED line into a dictionary."""
    fields = line.strip().split('\t')
    if 'KIR' not in fields[3]:
        return None
    attributes = fields[4].split(';')  # Assuming multiple attributes might be separated by ';'
    attributes_dict = {}
    for attr in attributes:
        key_value = attr.split('=')
        if len(key_value) == 2:
            attributes_dict[key_value[0]] = key_value[1]
    try:
        if '_e' in fields[3]:
            gene = fields[3].split('*')[0]
            allele = fields[3]
        elif '*' in fields[3]:
            gene, allele = fields[3].split('*')
            allele = allele[:3]
        elif '#' in fields[3]:
            gene = fields[3].split('#')[0]
            allele = fields[3][7:]
        else:
            raise ValueError(f"Error parsing gene and allele from {fields[3]}")
    except ValueError:
        print(f"Error parsing gene and allele from {fields[3]}")
        print(f"Line: {line}")
        raise

    return {
        'contig': fields[0],
        'source': 'IPD-KIR-V2.12.0',  # Placeholder source, update as needed
        'start': int(fields[1]),
        'end': int(fields[2]),
        'strand': fields[5],
        'attributes': {
            'template_allele': gene+ '*' + allele,
            'gene_name': gene
        }
    }

def load_and_convert_bed(file_path: str) -> t.List[t.Dict[str, t.Any]]:
    """Load a BED file and convert it into a list of dictionaries."""
    result = []
    with open(file_path, 'r') as file:
        for line in file:
            # Skip header lines and empty lines
            if line.startswith('track') or not line.strip():
                continue
            line_data = parse_bed_line(line)
            if line_data:
                result.append(line_data)

    return result


from typing import List
from collections import OrderedDict
import io

def load_csv_to_ordereddicts(file_path: str) -> List[OrderedDict]:
    """
    Loads data from a CSV file into a list of OrderedDicts.

    Args:
        file_path (str): The file path of the CSV to load.

    Returns:
        List[OrderedDict]: A list of OrderedDicts, each corresponding to a row in the CSV.
    """
    result = []  # List to hold the OrderedDicts
    with io.open(file_path, mode='r', encoding='utf-8') as file:
        # Read the first line to get the column names (headers)
        headers = [header.strip() for header in file.readline().split(',')]
        # Iterate through each remaining line in the file
        for line in file:
            # Split the line into values and strip whitespace
            values = [value.strip() for value in line.split(',')]
            # Combine the headers and values into an OrderedDict
            row_dict = OrderedDict(zip(headers, values))
            # Append the OrderedDict to the result list
            result.append(row_dict)
    return result


import msa_wrappers

from bakir import common
from bakir.common import Seq

from panel.pane import Bokeh


from importlib import reload

msa_aligner = msa_wrappers.ClustalWWrapper()

import yaml
from glob import glob
from bakir.alignment_variants import identify_closest_allele, call_variants, align_to_wildtype, identify_functional_variants
from IPython.display import Markdown
md = lambda x: display(Markdown(x))

def get_bakir_data(sample, haplo):
    with open(glob(f'HPRC-assemblies-annotations/{haplo}/{haplo}*.yaml')[0], 'r') as f:
        return yaml.load(f, Loader=yaml.UnsafeLoader)


def analyze_sample_gene_discordance(haplo, gene, db, alt_label='Immunanot', save_png=None):

    sample = haplo.split('.')[0]

    ka_geno = get_bakir_data(sample, haplo)
    
    ka_geno = group_annotations_by_gene(ka_geno)[gene]
    
    if alt_label == 'Immunanot':
        alt = group_annotations_by_gene_immunanot(extract_immunannot_gene_features(f'HPRC-Immunanot-annotations/{haplo}/{haplo}.gtf.gz'))[gene]
    elif alt_label == 'Skirt':
        alt = group_annotations_by_gene_immunanot(load_skirt_data(glob(f'HPRC-Skirt-annotations/{haplo}*/{haplo}*fa/{haplo}*allele.csv')[0]))[gene]
    else:
        raise ValueError("alt must be either 'immunanot' or 'skirt'")

    match_indexes, anno_matches, extra_kir, extra_imm = make_annotation_matches(ka_geno, alt)
    
    for c in extra_kir:
        md("#### Extra KA copy, showing:")
        analyze_extra_copy(ka_geno[c], db)

    for ka_anno, im_anno in anno_matches:
        ka_allele = ka_anno['closest allele'].split('*')[1][:3]
        im_allele = im_anno['template_allele'].split('*')[1][:3]

        if ka_allele != im_allele:
            closest_alleles, v, fv = analyze_discordant_alleles(ka_anno, im_anno, db=db, alt_label=alt_label, save_png=save_png)
            return closest_alleles, v, fv
    return None, None, None


import os
import contextlib

@contextlib.contextmanager
def suppress_output():
    original_stdout = os.dup(1)  # Duplicate the original stdout
    original_stderr = os.dup(2)  # Duplicate the original stderr
    try:
        with open(os.devnull, 'w') as devnull:
            os.dup2(devnull.fileno(), 1)  # Replace the current stdout with devnull
            os.dup2(devnull.fileno(), 2)  # Replace the current stderr with devnull
        yield
    finally:
        os.dup2(original_stdout, 1)  # Restore the original stdout
        os.dup2(original_stderr, 2)  # Restore the original stderr
        os.close(original_stdout)
        os.close(original_stderr)


def analyze_extra_copy(annotation, db):
    gene = annotation['gene']
    allele = annotation['closest allele'].split('*')[1]
    with suppress_output():
        alignment = msa_aligner.align([Seq('wild', db[gene].wildtype.seq), Seq('assembly', annotation['sequence']), Seq(allele, db[gene].alleles[allele].seq)])
    display(Bokeh(view_alignment(alignment, plot_width=1500)))
    
    
    align_records = {}
    for r in alignment:
        align_records[r.id] = r
        
    edit_distance = calc_ed(align_records['assembly'], align_records[allele])
    identity = 1 - (edit_distance / len(align_records['assembly']))
    
    if identity > 0.75:
        md(f"Copy valid: identity = {round(identity*100, 2)}%")
    else:
        md(f"Copy INVALID: identity = {round(identity*100, 2)}%")
                

def align_gene(gene, ka_allele, im_allele, assembly_seq, alt_assembly_seq=None, alt_label='Immunanot', save_png=None, display_alignment=True):
    short_labels = {"Immunanot": 'imm', 'Skirt': 'skirt'}
    seqs_to_align = [Seq('wild', gene.wildtype.seq), Seq('assembly', assembly_seq), 
                                       Seq(f'bk_allele-{ka_allele}', gene.alleles[ka_allele].seq),
                                       Seq(f'{short_labels[alt_label]}_allele-{im_allele}', gene.alleles[im_allele].seq)]
    if alt_assembly_seq:
        seqs_to_align.append(Seq(f'{short_labels[alt_label]}-assembly', alt_assembly_seq))
    with suppress_output():
        alignment = msa_aligner.align(seqs_to_align)
    if display_alignment:
        display(Bokeh(view_alignment(alignment, plot_width=1500, save_png=save_png)))
    return alignment

            
def calc_ed(seq1, seq2):
    ed = 0
    for i, j in zip(seq1.seq, seq2.seq):
        if i != j:
            ed += 1
    return ed


from Bio.Seq import Seq as BSeq

def analyze_discordant_alleles(ka_anno, im_anno, db, alt_label='Immunanot', save_png=None, display_table=True):
    gene = ka_anno['gene']
    
    ka_allele = ka_anno['closest allele']
    im_allele = im_anno['template_allele']
    
    if display_table:
        md(f"#### Discordant alleles: KA = {ka_allele.split('*')[1]}, {alt_label} = {im_allele.split('*')[1]}")

    
    # md(markdown_table(('source', 'allele', 'start', 'end'), []
    
    wildtype_sequence = db[gene].wildtype.seq

    poss_gene_seq, cigar_list = align_to_wildtype(ka_anno['sequence'], wildtype_sequence, seq_is_reversed=False)

    variants = call_variants(poss_gene_seq, cigar_list, wildtype_sequence)

    functional_variants = identify_functional_variants(variants, db[gene].wildtype)

    closest_alleles = pd.DataFrame(identify_closest_allele(variants, functional_variants, db[gene]))
    
    closest_alleles['common mut'] = closest_alleles['common mut'].apply(len)
    closest_alleles['new mut'] = closest_alleles['new mut'].apply(len)
    closest_alleles['missing mut'] = closest_alleles['missing mut'].apply(len)
    
    alt_assembly_seq = None
    if ka_anno['start'] != im_anno['start'] or ka_anno['end'] != im_anno['end']:
        whole_assembly = fetch_assembly_seq(ka_anno['sample'], ka_anno['haplotype'])
        alt_assembly_seq = str(whole_assembly[im_anno['contig']][im_anno['start']:im_anno['end']])
        if im_anno['strand'] == '-':
            alt_assembly_seq = str(BSeq(alt_assembly_seq).reverse_complement())
    alignment = align_gene(db[gene], ka_allele.split('*')[1], im_allele.split('*')[1], ka_anno['sequence'], alt_assembly_seq, alt_label=alt_label, save_png=save_png, display_alignment=display_table)
    
    align_records = {}
    for r in alignment:
        align_records[r.id] = r
    
    closest_alleles_filter = closest_alleles[closest_alleles['closest allele'].apply(lambda x: x in {ka_allele, im_allele})].copy()
    closest_alleles_filter['MSA edit distance'] = closest_alleles_filter['closest allele'].apply(lambda x: calc_ed(align_records[[y for y in align_records if x.split('*')[1] in y][0]], align_records['assembly']))
    
    if display_table:
        display(closest_alleles_filter)
    
    return closest_alleles_filter, variants, functional_variants

import pyfastx
def fetch_assembly_seq(sample, haplo):
    return pyfastx.Fasta(glob(f'HPRC-assemblies/{sample}/{sample}.{haplo}*.fa.gz')[0])


import panel as pn
import panel.widgets as pnw
pn.extension()
from bokeh.models import ColumnDataSource, Range1d, Rect, Text, FixedTicker, Span
from bokeh.plotting import figure
from bokeh.io import export_png  # Import the export_png function

def view_alignment(aln, fontsize="9pt", plot_width=800, save_png=None):
    """
    Generate and optionally save as PNG a Bokeh plot visualizing a sequence alignment.

    Args:
        aln (Iterable): An iterable of sequence records, each with `seq` and `id` attributes.
        fontsize (str, optional): Font size for sequence IDs. Defaults to "9pt".
        plot_width (int, optional): Width of the plot in pixels. Defaults to 800.
        save_png (str, optional): Path to save the plot as a PNG file. If None, the plot is not saved.

    Returns:
        Bokeh figure: The generated Bokeh plot object.
    """
    def get_colors(seqs):
        """Generate colors for bases in sequence."""
        text = [i for s in list(seqs) for i in s]
        clrs = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', '-': 'white', 'N': 'black'}
        colors = [clrs[i] for i in text]
        return colors

    # Prepare sequence and id lists from the alignment object
    seqs = [rec.seq for rec in aln]
    ids = [rec.id for rec in aln]
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)
    N = len(seqs[0])
    S = len(seqs)

    x = np.arange(1, N + 1)
    y = np.arange(0.5, S, 1)
    xx, yy = np.meshgrid(x, y)
    gx = xx.ravel()
    gy = yy.flatten()
    recty = gy

    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs) * 15 + 50
    x_range = Range1d(0, N + 1, bounds='auto')

    left_margin = max(len(id_) for id_ in ids) * 8

    p = figure(title=None, width=plot_width, height=plot_height,
               x_range=x_range, y_range=(0, S), tools="xpan,xwheel_zoom,reset,save",
               min_border=0, toolbar_location='below', min_border_left=left_margin)
    rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)

    tick_positions = [i + 0.5 for i in range(S)]
    p.yaxis.ticker = FixedTicker(ticks=tick_positions)
    p.yaxis.major_label_overrides = {i + 0.5: id_ for i, id_ in enumerate(ids)}
    p.yaxis.minor_tick_line_alpha = 0
    p.grid.visible = False
    p.yaxis.major_label_text_font_size = fontsize
    p.yaxis.axis_label_text_font_style = 'normal'

    for i in range(1, S):
        hline = Span(location=i, dimension='width', line_color='gray', line_width=1)
        p.add_layout(hline)

    # Save the plot as PNG if a path is provided
    print(f"Saving plot as PNG to {save_png}")
    if save_png:
        export_png(p, filename=save_png)

    return p

def markdown_table(header, data, transpose=False):
    '''Returns a string with table in markdown format
    Args
        header []           iterable of the header fields
        data []             iterable of data values (if 1 column) or of iterables containing column values (if > 1 column)
        transpose           transposes the table
    '''
    def table(header, data):
        result = []
        s = '|'+' | '.join(['{}'] * len(header))+'|'
        p = lambda x: s.format(*[str(i) for i in x])
        result.extend(list(map(p, (header if not isinstance(header[0], str) else [header]))))
        result.append('|'+('--|'*len(header)))
        if (not hasattr(data[0], '__iter__') or isinstance(data[0], str)):
            data = [[x] for x in data]
        result.extend(list(map(p, data)))
        return '\n'.join(result)

    if transpose:
        combined = [header] + data
        transposed = list(zip(*combined))
        return table(transposed[0], transposed[1:])
    
    return table(header, data)

def latex_table(header, data, transpose=False):
    def print_table(header, data):
        s = ' & '.join(['{: <20}'] * len(data[0])) + '\\\\'
        p = lambda x: s.format(*[str(i) for i in x])
        result = ['\\hline']
        result.extend(list(map(p, (header if not isinstance(header[0], str) else [header]))))
        result.append('\\hline')
        result.extend(list(map(p, data)))
        result.append('\\hline')
        print('\n'.join(result))
    
    if transpose:
        combined = [header] + data
        transposed = list(zip(*combined))
        print_table(transposed[0], transposed[1:])
    else:
        print_table(header, data)
