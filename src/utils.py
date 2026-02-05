import sys
from Bio import SeqIO
    
def clean_sequence(sequence):
    return str(sequence).strip().replace(" ", "").replace("\n", "").upper()

def validate_sequence(sequence):
    sequence = clean_sequence(sequence)
    IUPAC_letters = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'
    if all(aa in IUPAC_letters for aa in sequence):
        return True
    else:
        invalid_chr = set(sequence) - set(IUPAC_letters)
        print(f"Invalid characters in sequence: {', '.join(invalid_chr)}")
        return False

def validate_fasta(input_fasta, web=False):
    fh = None
    try:
        if web:
            from io import StringIO
            fh = StringIO(input_fasta)
        else:
            fh = open(input_fasta, "r")
            
        with fh:
            fasta_record = list(SeqIO.parse(fh, "fasta"))
            if len(fasta_record) == 0:
                return False
            for record in fasta_record:
                if not validate_sequence(record.seq):
                    return False
            return True
    except IOError:
        print(f"Error: The file could not be read.")
        return False
    except Exception as e:
        print(f"Error: {str(e)}, check if the file is in FASTA format.")
        return False

def format_result(inputs, results, format="txt"):
    if format == "txt":
        result = f"pEpitope calculator result:\n" + "-"* 60 + "\n" + \
        f"Dominat strain: {inputs['name']}\n" + \
        f"Vaccine strain: {inputs['vaccine_name']}\n" + \
        f"Subtype: {inputs['subtype']}\n" + \
        f"Dataset/Mode: {inputs['dataset']}/{inputs['mode']}\n\nResults\n" + \
        f"- pEpitope: {results['pepitope_score']}\n" + \
        f"- Dominant epitope: {results['dominant_epitope']}\n" + \
        f"- Predicted VE: {results['ve']}\n" + \
        f"- STD Error: {results['ve_std_err']}\n\n  Substitutions\n" + \
        "\n".join(f"  - Epitope {i}: {', '.join(results['substitutions'][i])}" for i in ["A", "B", "C", "D", "E"])
        return result
    elif format == "json":
        import json
        result = {
            "subtype": inputs["subtype"],
            "dataset": inputs["dataset"],
            "mode": inputs["mode"],
            "dominant_strain": {
                "name": inputs["name"]
            },
            "vaccine_strain": {
                "name": inputs["vaccine_name"]
            },
            "results": {
                "pepitope_score": results["pepitope_score"],
                "dominant_epitope": results["dominant_epitope"],
                "ve": results["ve"],
                "ve_std_err": results["ve_std_err"]
            },
            "substitutions": {
                "A": results["substitutions"]["A"],
                "B": results["substitutions"]["B"],
                "C": results["substitutions"]["C"],
                "D": results["substitutions"]["D"],
                "E": results["substitutions"]["E"]
            }
        }
        return json.dumps(result, indent=2)
    else:
        raise ValueError(f"Invalid format: {format}")