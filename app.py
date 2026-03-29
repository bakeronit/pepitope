import asyncio
import json
from pyscript import when, document
from src.core import generate_results
from src.models import get_model

EXAMPLE_DOMINANT = """>A/Switzerland/9715293/2013
MKAKLLVLLYAFVATDADTQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWAGVTQNGTSSSCIRGSNSSFFSRLNWLTHLNSKYPALNVTMPNNEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPERQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGRGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAMSCFLLCVALLGFIMWACQKGSLQCRICI"""

EXAMPLE_VACCINE = """>A/Hong Kong/4801/2014 (egg-adapted)
MKTIIALSYILCLVFAQKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSSSSFFSRLNWLTHLNYTYPALNVTMPNNEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKHSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGRGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHNVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"""

btn = document.querySelector("#calculate")
btn_label = document.querySelector("#btn-label")
btn_spinner = document.querySelector("#btn-spinner")

btn.disabled = False
btn_label.textContent = "Calculate"
btn_spinner.style.display = "none"

document.querySelector("#dominant").value = EXAMPLE_DOMINANT
document.querySelector("#vaccine").value = EXAMPLE_VACCINE

output = document.querySelector("#output")


def render_result(data):
    r = data["results"]
    subs = data["substitutions"]

    ve = f"{float(r['ve'].replace(' %', '')):.2f} %"
    ve_std_err = f"{float(r['ve_std_err'].replace(' %', '')):.2f} %"
    pepitope_score = f"{float(r['pepitope_score']):.4f}"

    sub_rows = ""
    for epitope in ["A", "B", "C", "D", "E"]:
        items = subs.get(epitope, [])
        display = ", ".join(items) if items else "\u2014"
        sub_rows += f"<tr><td>Epitope {epitope}</td><td>{display}</td></tr>"

    return f"""
    <div class="result-header">
        <h3>Results</h3>
        <span class="result-badge">{data['subtype']} &middot; {data['dataset']} &middot; {data['mode']}</span>
    </div>
    <div class="result-strains">
        <div><span class="strain-label">Dominant</span>{data['dominant_strain']['name']}</div>
        <div><span class="strain-label">Vaccine</span>{data['vaccine_strain']['name']}</div>
    </div>
    <div class="result-metrics">
        <div class="metric">
            <div class="metric-value">{ve}</div>
            <div class="metric-label">Predicted VE</div>
        </div>
        <div class="metric">
            <div class="metric-value">{pepitope_score}</div>
            <div class="metric-label">pEpitope Score</div>
        </div>
        <div class="metric">
            <div class="metric-value">{r['dominant_epitope']}</div>
            <div class="metric-label">Dominant Epitope</div>
        </div>
        <div class="metric">
            <div class="metric-value">{ve_std_err}</div>
            <div class="metric-label">Std Error</div>
        </div>
    </div>
    <div class="result-substitutions">
        <h4>Substitutions</h4>
        <table>
            <tbody>{sub_rows}</tbody>
        </table>
    </div>
    """


@when("click", "#calculate")
async def calculate(event) -> None:
    dominant_sequence = document.querySelector("#dominant").value
    vaccine_sequence = document.querySelector("#vaccine").value
    vaccine_type = document.querySelector("[name='vaccine_type']:checked").value
    dataset = document.querySelector("[name='dataset']:checked").value
    mode = document.querySelector("[name='mode']:checked").value

    output.innerHTML = ""
    output.className = ""
    output.style.display = "none"

    btn.disabled = True
    btn_label.textContent = "Calculating…"
    btn_spinner.style.display = "inline-block"
    await asyncio.sleep(0.01)

    try:
        model = get_model(vaccine_type)
        if not dominant_sequence:
            raise ValueError("Please enter a dominant sequence.")
        if not vaccine_sequence:
            raise ValueError("Please enter a vaccine sequence.")
        result_json = generate_results(
            dominant_sequence, vaccine_sequence, model, dataset, mode,
            outfmt="json", from_string=True
        )
        records = json.loads(result_json)
        if isinstance(records, dict):
            records = [records]
        output.className = "result-success"
        output.innerHTML = "".join(render_result(d) for d in records)
        output.style.display = "block"
    except Exception as e:
        output.className = "result-error"
        output.innerHTML = f'<p class="error-msg">{str(e)}</p>'
        output.style.display = "block"
    finally:
        btn.disabled = False
        btn_label.textContent = "Calculate"
        btn_spinner.style.display = "none"
