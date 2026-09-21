#!/usr/bin/env python3
"""Generate optimized repeat-unit SMILES from prepared ICL prompts."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd


def render_prompt(record: dict[str, object]) -> str:
    examples = "\n".join(
        f"Original SMILES: {item['SMILES_original']}\nOptimized SMILES: {item['SMILES_optimized']}"
        for item in record["examples"]
    )
    question = record["question"]
    return (
        f"{record['context']}\n\n{examples}\n\n{question['instruction']}\n"
        f"Original SMILES: {question['SMILES_test']}\n"
        f"Expected output: {record['expected_output_format']}"
    )


def openai_generator(model: str):
    if not os.getenv("OPENAI_API_KEY"):
        raise RuntimeError("Set OPENAI_API_KEY before using an OpenAI model")
    from openai import OpenAI

    client = OpenAI()

    def generate(prompt: str) -> str:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": "You are a polymer chemistry expert."},
                {"role": "user", "content": prompt},
            ],
            temperature=1.0,
            max_tokens=512,
        )
        return response.choices[0].message.content.strip()

    return generate


def hf_generator(model: str):
    from transformers import pipeline

    pipe = pipeline("text-generation", model=model, device_map="auto")

    def generate(prompt: str) -> str:
        result = pipe(prompt, max_new_tokens=512, do_sample=False, return_full_text=False)
        return result[0]["generated_text"].strip()

    return generate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prompts", type=Path, required=True)
    parser.add_argument("--provider", choices=["openai", "huggingface"], required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = json.loads(args.prompts.read_text(encoding="utf-8"))
    generate = openai_generator(args.model) if args.provider == "openai" else hf_generator(args.model)
    rows = []
    for pid, record in records.items():
        response = generate(render_prompt(record))
        optimized = response.split("Optimized SMILES:", 1)[-1].strip()
        rows.append(
            {
                "PID": pid,
                "Original_SMILES": record["question"]["SMILES_test"],
                "Optimized_SMILES": optimized,
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(f"Wrote {len(rows)} generations to {args.output}")


if __name__ == "__main__":
    main()

