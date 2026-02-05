import sys
from Bio.Align import PairwiseAligner
from Bio import SeqIO

from .calculation import calculate_pepitope_score, calculate_vaccine_efficacy, extract_substitutions
from .utils import validate_fasta, format_result

aligner = PairwiseAligner()
aligner.mode = "local"

def get_optimal_alignment(seq, ref_seq):
    alignment = aligner.align(ref_seq, seq)[0]
    return alignment

def get_epitope_aa(
    alignment,
    epitope_positions: list[int],
    ref_offset: int 
    ) -> list[str]:
    target_positions = [ pos + ref_offset - 1 for pos in epitope_positions ]
    epitope_aa  = []
    n_gap = 0
    for align_pos, letter in enumerate(alignment[0]):
        if letter == '-':
            n_gap += 1
            target_positions = [ pos + 1 if pos >= align_pos else pos for pos in target_positions ]
            if align_pos >= min(target_positions) and align_pos <= max(target_positions):
                raise ValueError("Error: gap is in the epitope region")
        elif align_pos in target_positions:
            epitope_aa.append(f"{letter}{align_pos - n_gap - ref_offset + 1}{alignment[1][align_pos]}")
    
    return epitope_aa


def generate_results(dominant_fasta, vaccine_fasta, model, dataset, mode, outfmt="txt", web = False):
    if not validate_fasta(dominant_fasta, web):
        raise ValueError("Invalid dominant FASTA file.")
    
    if not validate_fasta(vaccine_fasta, web):
        raise ValueError("Invalid vaccine FASTA file.")

    dataset = dataset.upper()
    mode = mode.capitalize()
    
    if (dataset, mode) not in model.ve_equation_params:
         raise ValueError(f"Invalid dataset/mode combination: {dataset}/{mode}. Available: {list(model.ve_equation_params.keys())}")

    if web:
        from io import StringIO
        vaccine_fasta = StringIO(vaccine_fasta)
    else:
        vaccine_fasta = open(vaccine_fasta, "r")

    with vaccine_fasta as fh:
        vaccine_strain = SeqIO.read(fh, "fasta")

    vaccine_alignment = get_optimal_alignment(vaccine_strain, model.sequence)
    vaccine_epitopes = {
        site : get_epitope_aa(vaccine_alignment, model.epitope_positions[site], model.ha1_region[0])
        for site in model.epitope_positions.keys()
    }

    ve_params = model.ve_equation_params[(dataset, mode)]
    if web:
        dominant_fh = StringIO(dominant_fasta)
    else:
        dominant_fh = open(dominant_fasta, "r")

    with dominant_fh:
        for record in SeqIO.parse(dominant_fh, "fasta"):
            alignment = get_optimal_alignment(record, model.sequence)
            epitopes = {
                site : get_epitope_aa(alignment, model.epitope_positions[site], model.ha1_region[0])
                for site in model.epitope_positions.keys()
            }
            pepitope_scores = {
                site : calculate_pepitope_score(epitopes[site], vaccine_epitopes[site])
                for site in model.epitope_positions.keys()
            }
            epitope_substitutions = {
                site : extract_substitutions(epitopes[site], vaccine_epitopes[site])
                for site in model.epitope_positions.keys()
            }

            dominant_epitope = max(pepitope_scores, key=pepitope_scores.get)
            vaccine_efficacy, vaccine_efficacy_err = calculate_vaccine_efficacy(pepitope_scores[dominant_epitope], ve_params)

            inputs = {
                "name" : record.description,
                "vaccine_name": vaccine_strain.description,
                "subtype" : model.subtype,
                "dataset" : dataset,
                "mode" : mode,
            }

            results = {
                "pepitope_score" : f"{pepitope_scores[dominant_epitope]:.6f}",
                "dominant_epitope" : dominant_epitope,
                "ve" : f"{vaccine_efficacy * 100:.4f} %",
                "ve_std_err" : f"{vaccine_efficacy_err * 100:.4f} %",
                "substitutions" : epitope_substitutions
            }

            return format_result(inputs, results, outfmt)