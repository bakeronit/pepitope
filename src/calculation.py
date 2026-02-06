from math import sqrt

def calculate_pepitope_score(epitope: list, vaccine_epitope: list) -> float:
    if len(epitope) != len(vaccine_epitope):
        raise ValueError("Error: epitope and vaccine epitope must have the same length")
    return sum(1 for a, b in zip(epitope, vaccine_epitope) if a != b) / len(epitope)

def calculate_vaccine_efficacy(pepitope_score: float, params: dict) -> list[float, float]:
    m = params['m']
    m_err = params['m_err']
    b = params['b']
    b_err = params['b_err']
    return [m * pepitope_score + b, sqrt(m_err ** 2 * pepitope_score ** 2 + b_err ** 2)]

def extract_substitutions(epitope: list, vaccine_epitope: list) -> list[str]:
    return [f"{b[-1]}{a[1:]}" for a, b in zip(epitope, vaccine_epitope) if a != b]