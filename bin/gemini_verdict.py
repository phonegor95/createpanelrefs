#!/usr/bin/env python3
"""
Independent LLM second-opinion on a per-sample HBA alpha-thalassemia call.

Sends ONE sample's HBA_CLASSIFIER YAML (copy numbers + the full evidence block:
depth, allelic-LOH, junction split-reads/soft-clips) to the Gemini API and asks
for a terse, structured VERDICT on whether the GATK-based call is supported by
its own evidence -- a downstream adversarial check, NOT ground truth.

Design notes
------------
* The API key is read from the environment (GEMINI_API_KEY), never from argv,
  so it does not land in Nextflow's .command.sh. Under Nextflow use the
  `secret 'GEMINI_API_KEY'` directive; standalone, export it first.
* Network egress is optional via --proxy (the box reaches the API through a
  local proxy). stdlib only (urllib) so it runs in the GATK container as-is.
* Advisory module: on any API/network error it still writes a verdict file
  (status: ERROR) and exits 0, so it never breaks the pipeline. The call itself
  is unchanged -- this only annotates.
"""
import argparse
import json
import os
import sys
import urllib.request
import urllib.error

PROMPT_HEADER = """You are a clinical genomics reviewer. Below is ONE sample's automated \
alpha-globin (HBA) classifier output (YAML). The caller fuses three orthogonal axes over a \
recal.cram: region DEPTH, allelic loss-of-heterozygosity, and junction split-reads/soft-clips. \
Rule: a deletion (CN=1) is called only if depth and junction corroborate; a lone junction with \
copy-neutral depth is flagged "structural variant, manual review"; the allelic axis is \
supporting-only. Note -a4.2 is depth- and junction-silent in the segdup (a known blind spot).

Give a terse VERDICT in this exact format:
VERDICT: <AGREE | AGREE-WITH-CAVEATS | DISAGREE>
CONFIDENCE: <high | medium | low>
REASONING: <2-4 sentences tying the call to the evidence numbers shown>
WATCH: <one concrete thing to double-check, or "none">

Sample YAML:
"""


def build_prompt(yaml_text):
    return PROMPT_HEADER + "\n" + yaml_text.strip() + "\n"


def call_gemini(prompt, model, api_key, proxy, timeout=120):
    url = (f"https://generativelanguage.googleapis.com/v1beta/models/"
           f"{model}:generateContent?key={api_key}")
    body = json.dumps({
        "contents": [{"parts": [{"text": prompt}]}],
        # gemini-3.5-flash is a thinking model: thoughtsTokenCount counts against
        # maxOutputTokens, so the cap must cover internal reasoning + the answer.
        "generationConfig": {"temperature": 0.2, "maxOutputTokens": 2048},
    }).encode()
    req = urllib.request.Request(url, data=body,
                                 headers={"Content-Type": "application/json"})
    if proxy:
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({"https": proxy, "http": proxy}))
    else:
        opener = urllib.request.build_opener()
    r = opener.open(req, timeout=timeout)
    d = json.load(r)
    text = d["candidates"][0]["content"]["parts"][0]["text"]
    um = d.get("usageMetadata", {})
    return text, um


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--yaml", required=True, help="per-sample HBA.yaml")
    ap.add_argument("--model", default="gemini-3.5-flash")
    ap.add_argument("--proxy", default="", help="http(s) proxy URL, optional")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    yaml_text = open(a.yaml).read()
    api_key = os.environ.get("GEMINI_API_KEY", "").strip()

    header = (f"# HBA independent verdict  sample={a.sample}\n"
              f"# reviewer={a.model} (downstream second-opinion, NOT ground truth)\n")

    if not api_key:
        open(a.out, "w").write(header + "status: ERROR\nVERDICT: NA\n"
                               "REASONING: GEMINI_API_KEY not set in environment.\n")
        sys.stderr.write("WARN: GEMINI_API_KEY missing; wrote ERROR verdict.\n")
        return

    try:
        text, um = call_gemini(build_prompt(yaml_text), a.model, api_key, a.proxy)
        usage = (f"# tokens: prompt={um.get('promptTokenCount')} "
                 f"out={um.get('candidatesTokenCount')}\n")
        open(a.out, "w").write(header + "status: OK\n" + usage + "\n" + text.strip() + "\n")
        print(f"{a.sample}: gemini verdict OK")
    except urllib.error.HTTPError as e:
        msg = e.read().decode()[:300]
        open(a.out, "w").write(header + f"status: ERROR\nVERDICT: NA\n"
                               f"REASONING: HTTP {e.code} {msg}\n")
        sys.stderr.write(f"WARN: gemini HTTP {e.code}; wrote ERROR verdict.\n")
    except Exception as e:  # network/timeout/parse -- stay advisory, never fail the run
        open(a.out, "w").write(header + f"status: ERROR\nVERDICT: NA\n"
                               f"REASONING: {type(e).__name__}: {e}\n")
        sys.stderr.write(f"WARN: gemini call failed ({type(e).__name__}); wrote ERROR verdict.\n")


if __name__ == "__main__":
    main()
