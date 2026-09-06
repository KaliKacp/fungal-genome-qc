"""
Renderuje finalny raport HTML na podstawie results/summary.json
i szablonu Jinja2.
"""
import json
import yaml
from jinja2 import Environment, FileSystemLoader

def load_config(path: str = "config/config.yaml") -> dict:
    with open(path) as f:
        return yaml.safe_load(f)

def load_summary(path: str = "results/summary.json") -> dict:
    with open(path) as f:
        return json.load(f)

def render(config: dict, summary: dict, template_dir: str = "workflow/templates") -> str:
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("report_template.html")

    return template.render(
        organism=config["organism"],
        accession=config["assembly_accession"],
        seqkit=summary["seqkit"],
        quast=summary["quast"],
        busco=summary.get("busco", {"lineage": config["busco_lineage"]}),
    )

if __name__ == "__main__":
    config = load_config()
    summary = load_summary()
    html = render(config, summary)

    with open("results/report.html", "w") as f:
        f.write(html)

    print("Raport zapisany: results/report.html")