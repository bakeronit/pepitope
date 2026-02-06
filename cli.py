import typer
from src import core

def main(dominant_fasta: str, vaccine_fasta: str, vaccine_type: str = "H3N2", dataset: str = "CDC", mode: str = "Efficacy", outfmt: str = "txt"):
    if vaccine_type.upper() not in ["H3N2"]:
        raise ValueError("Invalid vaccine type")
    else:
        from src.model import H3N2
        model = H3N2
    dataset = dataset.upper()
    if dataset not in ["CDC", "NH"]:
        raise ValueError("Invalid dataset")
    
    mode = mode.capitalize()
    if mode not in ["Efficacy", "Effectiveness"]:
        raise ValueError("Invalid mode")
    
    try:
        result = core.generate_results(dominant_fasta, vaccine_fasta, model, dataset, mode, outfmt)
        typer.echo(result)
    except Exception as e:
        typer.secho(f"Error: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

if __name__ == "__main__":
    typer.run(main)