import typer
from src import core
from src.models import get_model

def main(dominant_fasta: str, vaccine_fasta: str, vaccine_type: str = "H3N2", dataset: str = "CDC", mode: str = "Efficacy", outfmt: str = "txt"):
    model = get_model(vaccine_type)
    dataset = dataset.upper()
    mode = mode.capitalize()
    
    try:
        result = core.generate_results(dominant_fasta, vaccine_fasta, model, dataset, mode, outfmt)
        typer.echo(result)
    except Exception as e:
        typer.secho(f"Error: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

if __name__ == "__main__":
    typer.run(main)