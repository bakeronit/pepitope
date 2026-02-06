import asyncio
from pyscript import when, document
from src.core import generate_results
from src.model import H3N2

btn = document.querySelector("#calculate")
btn.disabled = False

output = document.querySelector("#output")

@when("click", "#calculate")
async def calculate(event) -> None:
    dominant_sequence = document.querySelector("#dominant").value
    vaccine_sequence = document.querySelector("#vaccine").value
    vaccine_type = document.querySelector("[name='vaccine_type']:checked").value
    dataset = document.querySelector("[name='dataset']:checked").value
    mode = document.querySelector("[name='mode']:checked").value

    output.textContent = ""
    output.style.display = "none"
    await asyncio.sleep(0.01)

    if vaccine_type == "h3n2":
        model = H3N2
    else:
        ...  # future work to implement other vaccine types
    
    try:
        if not dominant_sequence:
            raise ValueError("Please enter a dominant sequence.")
        if not vaccine_sequence:
            raise ValueError("Please enter a vaccine sequence.")
        results = generate_results(dominant_sequence, vaccine_sequence, model, dataset, mode, web=True)
        output.style.display = "block"
        output.textContent = results
    except Exception as e:
        output.style.display = "block"
        output.textContent = f"Error: {str(e)}"
    finally:
        btn.disabled = False