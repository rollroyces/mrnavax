"""LLM-calling wrapper.

Two backends, selected at runtime:

1. ``hermes``  — when running inside the Hermes desktop app the agent itself can
   simply invoke this module's ``llm_complete`` from its own context. The
   ``hermes`` backend just records the call (no nested LLM in the demo) and
   returns a deterministic mock so the tool still runs end-to-end.

2. ``openai``  — when the ``OPENAI_API_KEY`` env var is set, calls go to the
   OpenAI Chat Completions API at ``gpt-4o-mini`` (cheap + fast, fine for tool
   use). Any OpenAI-compatible endpoint can be selected via ``OPENAI_BASE_URL``.

3. ``mock``    — deterministic stub for tests and offline runs. Returns the
   prompt's last user message verbatim, prefixed with ``"[mock] "``.

The CLI lets the caller pick a backend with ``--backend {mock,openai,hermes}``.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import Any

BACKEND = os.environ.get("MRNA_AI_LLM_BACKEND", "auto").lower()


def _detect_backend() -> str:
    if BACKEND in {"mock", "openai", "hermes"}:
        return BACKEND
    # auto
    if os.environ.get("OPENAI_API_KEY"):
        return "openai"
    return "mock"


def llm_complete(
    prompt: str,
    *,
    system: str = "",
    json_mode: bool = False,
    temperature: float = 0.2,
    max_tokens: int = 800,
    backend: str | None = None,
) -> str:
    """Return a single LLM completion for ``prompt``.

    Backend selection: explicit ``backend`` arg > ``MRNA_AI_LLM_BACKEND`` env
    > auto-detect.
    """
    backend = backend or _detect_backend()

    if backend == "mock":
        # Deterministic offline stub: extract a structured guess from the prompt.
        return _mock_complete(prompt, system=system, json_mode=json_mode)

    if backend == "openai":
        return _openai_complete(
            prompt,
            system=system,
            json_mode=json_mode,
            temperature=temperature,
            max_tokens=max_tokens,
        )

    if backend == "hermes":
        # Inside the Hermes app, the assistant already IS the LLM. We surface
        # a clear error so the operator knows to switch backends.
        raise RuntimeError(
            "backend='hermes' is reserved for the assistant's own context. "
            "Use backend='openai' with OPENAI_API_KEY, or backend='mock' for "
            "offline runs."
        )

    raise ValueError(f"unknown backend: {backend!r}")


# ---------- mock backend ----------------------------------------------------


def _mock_complete(prompt: str, *, system: str = "", json_mode: bool = False) -> str:
    """Deterministic stub. Parses the prompt for obvious structured fields."""
    text = prompt.lower()
    if json_mode or "return json" in text or "json object" in text:
        # Heuristic: extract the first peptide or HLA mentioned.
        import re

        hla = re.search(r"hla[-_][a-z0-9*:.]+", prompt, re.I)
        pep = re.search(r"peptide\s*[:=]\s*([A-Z]{8,11})", prompt)
        var = re.search(r"variant\s*[:=]\s*([A-Z]\d+[A-Z])", prompt)
        if pep and hla:
            guess = {
                "peptide": pep.group(1),
                "hla": hla.group(0),
                "binding_affinity_nM": 250.0 if "A*02" in hla.group(0) else 800.0,
                "binder": True if "A*02" in hla.group(0) else False,
                "immunogenicity_score": 0.55,
                "rationale": "[mock] HLA-A*02:01 favors hydrophobic anchors at P2/P9.",
            }
            return json.dumps(guess)
        if var:
            return json.dumps(
                {
                    "variant": var.group(1),
                    "immunogenic": True,
                    "rationale": "[mock] missense variant in known tumor-suppressor locus.",
                }
            )
        # Trial-matching prompt: extract bullet-pointed criteria and judge each
        # by simple keyword overlap with the patient summary.
        # Stop parsing at the JSON shape section that follows.
        end_markers = ["return only a json object", "return json object", "json object with this"]
        # Better: split at the "Exclusion criteria" marker.
        if "exclusion criteria" in text:
            inc_block = prompt.split("Inclusion criteria")[1].split("Exclusion criteria")[0]
            exc_block = prompt.split("Exclusion criteria")[1]
            # Truncate exc_block at any "return json" or definition list
            for marker in end_markers:
                if marker in exc_block.lower():
                    exc_block = exc_block.lower().split(marker)[0]
        else:
            inc_block, exc_block = prompt, ""
        inc_items = re.findall(r"-\s*(.+?)(?:\n|$)", inc_block)
        exc_items = re.findall(r"-\s*(.+?)(?:\n|$)", exc_block)
        # Patient summary block (case-insensitive marker)
        if "patient summary" in text:
            idx = prompt.lower().find("patient summary")
            tail = prompt[idx + len("patient summary:") :]
            trial_idx = tail.lower().find("trial:")
            pat_block = tail[:trial_idx] if trial_idx >= 0 else tail
        else:
            pat_block = prompt
        pat_lower = pat_block.lower()

        def judge(criterion: str, is_exclusion: bool = False) -> str:
            crit_tokens = [t for t in re.split(r"[^a-z0-9]+", criterion.lower()) if len(t) > 3]
            if not crit_tokens:
                return "uncertain"
            # Negation in the patient text
            negations = ["no ", "not ", "denies ", "without "]
            pat_negated = any(neg in pat_lower for neg in negations)
            # Count hits: each token must appear, but for exclusions we
            # also check for negation in the patient text.
            hits = sum(1 for t in crit_tokens if t in pat_lower)
            # For exclusion criteria, a negation in the patient text
            # about a relevant concept inverts the meaning (e.g., "no
            # prior therapy" + criterion "Prior therapy" → "unmet").
            if is_exclusion and pat_negated:
                # Only flip if the negation is "close" to the criterion
                # concept (heuristic: any token of the criterion appears
                # after the negation).
                for neg in negations:
                    idx = pat_lower.find(neg)
                    while idx >= 0:
                        rest = pat_lower[idx + len(neg) : idx + len(neg) + 80]
                        if any(t in rest for t in crit_tokens):
                            return "unmet"
                        idx = pat_lower.find(neg, idx + 1)
            if hits >= max(1, len(crit_tokens) // 2):
                return "met"
            if is_exclusion and hits == 0:
                return "unmet"
            return "uncertain"

        return json.dumps(
            {
                "inclusion": [
                    {
                        "criterion": c,
                        "verdict": judge(c, is_exclusion=False),
                        "evidence": "[mock] keyword overlap",
                    }
                    for c in inc_items
                ],
                "exclusion": [
                    {
                        "criterion": c,
                        "verdict": judge(c, is_exclusion=True),
                        "evidence": "[mock] keyword overlap",
                    }
                    for c in exc_items
                ],
            }
        )
    return "[mock] " + prompt.strip().splitlines()[-1][:200]


# ---------- openai backend --------------------------------------------------


def _openai_complete(
    prompt: str,
    *,
    system: str,
    json_mode: bool,
    temperature: float,
    max_tokens: int,
) -> str:
    api_key = os.environ["OPENAI_API_KEY"]
    base = os.environ.get("OPENAI_BASE_URL", "https://api.openai.com/v1")
    model = os.environ.get("OPENAI_MODEL", "gpt-4o-mini")

    body: dict[str, Any] = {
        "model": model,
        "messages": ([{"role": "system", "content": system}] if system else [])
        + [{"role": "user", "content": prompt}],
        "temperature": temperature,
        "max_tokens": max_tokens,
    }
    if json_mode:
        body["response_format"] = {"type": "json_object"}

    req = urllib.request.Request(
        f"{base}/chat/completions",
        data=json.dumps(body).encode(),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"OpenAI HTTP {e.code}: {e.read().decode()[:500]}") from e

    try:
        return data["choices"][0]["message"]["content"]
    except (KeyError, IndexError) as e:
        raise RuntimeError(f"unexpected OpenAI response: {data}") from e


def llm_json(prompt: str, **kwargs: Any) -> dict[str, Any]:  # noqa: ANN401
    """Convenience: call ``llm_complete`` with json_mode=True and parse."""
    raw = llm_complete(prompt, json_mode=True, **kwargs)
    # Strip code fences if any
    raw = raw.strip()
    if raw.startswith("```"):
        raw = raw.split("\n", 1)[1].rsplit("```", 1)[0].strip()
    try:
        return json.loads(raw)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"LLM did not return valid JSON: {raw[:300]!r}") from e
