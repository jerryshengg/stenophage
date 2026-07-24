#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd


INPUT_FILE = Path("cog_annotations.csv")
OUTPUT_FILE = Path("cog_annotations_new.csv")


# COG category X:
# Mobilome, prophages, and transposons.
X_KEYWORDS = [
    r"\bphage\b",
    r"\bprophage\b",
    r"\bviral\b",
    r"\bvirion\b",
    r"\bcapsid\b",
    r"\btail protein\b",
    r"\btail fiber\b",
    r"\btail fibre\b",
    r"\btail sheath\b",
    r"\btail tube\b",
    r"\bbaseplate\b",
    r"\bportal protein\b",
    r"\bterminase\b",
    r"\bholin\b",
    r"\bendolysin\b",
    r"\blysin\b",
    r"\bintegrase\b",
    r"\btransposase\b",
    r"\btransposon\b",
    r"\bretrotransposon\b",
    r"\binsertion sequence\b",
    r"\bmobile element\b",
    r"\bmobilome\b",
]


def clean_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Remove whitespace and leading comment symbols from column names."""
    df = df.copy()
    df.columns = [
        str(column).strip().lstrip("#").strip()
        for column in df.columns
    ]
    return df


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Find a column without requiring exact capitalization."""
    normalized = {
        str(column).strip().lower().replace(" ", "_"): column
        for column in df.columns
    }

    for candidate in candidates:
        key = candidate.strip().lower().replace(" ", "_")
        if key in normalized:
            return normalized[key]

    return None


def clean_text(value: object) -> str:
    """Convert missing annotation values to an empty string."""
    if pd.isna(value):
        return ""

    text = str(value).strip()

    if text.lower() in {"", "-", "nan", "none", "na", "n/a"}:
        return ""

    return text


def combine_annotations(row: pd.Series, annotation_columns: list[str]) -> str:
    """
    Combine available annotation fields while avoiding duplicates.

    Example:
        Preferred_name: terL
        Description: Large terminase subunit

    Output:
        terL | Large terminase subunit
    """
    annotations = []

    for column in annotation_columns:
        text = clean_text(row.get(column, ""))

        if text and text not in annotations:
            annotations.append(text)

    return " | ".join(annotations)


def normalize_cog_category(value: object) -> str:
    """
    Retain only valid single-letter COG category characters.

    Examples:
        "E,G"  -> "EG"
        "MNO"  -> "MNO"
        "S"    -> "S"
        "-"    -> ""
    """
    text = clean_text(value).upper()

    if not text:
        return ""

    letters = re.findall(r"[A-Z]", text)

    # Preserve order while removing duplicate letters.
    return "".join(dict.fromkeys(letters))


def annotation_matches_x(annotation: str) -> bool:
    """Return True when an annotation supports COG category X."""
    return any(
        re.search(pattern, annotation, flags=re.IGNORECASE)
        for pattern in X_KEYWORDS
    )


def assign_x_when_supported(row: pd.Series) -> str:
    """
    Assign X only when no informative COG category exists and the annotation
    indicates a mobilome-, prophage-, transposon-, or virion-related function.
    """
    category = normalize_cog_category(row["COG_category_original"])
    annotation = clean_text(row["combined_annotation"])

    # Keep existing informative eggNOG categories unchanged.
    if category and category != "S":
        return category

    if annotation_matches_x(annotation):
        return "X"

    # Preserve S when eggNOG specifically assigned unknown function.
    if category == "S":
        return "S"

    return "-"


def split_cog_categories(df: pd.DataFrame) -> pd.DataFrame:
    """
    Expand multi-letter COG categories into one row per category.

    Example:
        gene_1, EG

    becomes:
        gene_1, E
        gene_1, G

    Unassigned rows remain as "-".
    """
    expanded_rows = []

    for _, row in df.iterrows():
        category = clean_text(row["COG_category"])

        if not category:
            category = "-"

        if category == "-":
            categories = ["-"]
        else:
            categories = list(dict.fromkeys(re.findall(r"[A-Z]", category)))

            if not categories:
                categories = ["-"]

        for single_category in categories:
            new_row = row.copy()
            new_row["COG_category"] = single_category
            expanded_rows.append(new_row)

    return pd.DataFrame(expanded_rows)


def main() -> None:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Could not find {INPUT_FILE.resolve()}\n"
            "Place the script in the same folder as cog_annotations.csv."
        )

    try:
        df = pd.read_csv(INPUT_FILE, dtype=str)
    except UnicodeDecodeError:
        df = pd.read_csv(INPUT_FILE, dtype=str, encoding="latin-1")
    except pd.errors.ParserError:
        # Fallback in case the file is actually tab-delimited despite .csv.
        df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str)

    df = clean_column_names(df)

    cog_column = find_column(
        df,
        [
            "COG_category",
            "COG category",
            "cog_category",
            "COG",
        ],
    )

    if cog_column is None:
        raise ValueError(
            "No COG category column was found.\n"
            f"Available columns: {list(df.columns)}"
        )

    preferred_name_column = find_column(
        df,
        [
            "Preferred_name",
            "preferred name",
            "preferred_name",
            "gene",
            "gene_name",
        ],
    )

    description_column = find_column(
        df,
        [
            "Description",
            "description",
            "annotation",
            "function",
            "product",
        ],
    )

    annotation_columns = [
        column
        for column in [preferred_name_column, description_column]
        if column is not None
    ]

    # Rename the detected field to a consistent name.
    if cog_column != "COG_category":
        df = df.rename(columns={cog_column: "COG_category"})

    df.insert(
        df.columns.get_loc("COG_category"),
        "COG_category_original",
        df["COG_category"].fillna("-"),
    )

    if annotation_columns:
        df["combined_annotation"] = df.apply(
            combine_annotations,
            axis=1,
            annotation_columns=annotation_columns,
        )
    else:
        df["combined_annotation"] = ""

    # Standardize the categories and assign X where biologically supported.
    df["COG_category"] = df.apply(assign_x_when_supported, axis=1)

    # Record whether X was newly assigned by this script.
    df["X_added"] = (
        df["COG_category"].eq("X")
        & ~df["COG_category_original"]
        .fillna("")
        .str.upper()
        .str.contains("X", regex=False)
    )

    # Split categories such as EG or MNO into separate rows.
    output_df = split_cog_categories(df)

    output_df.to_csv(OUTPUT_FILE, index=False)

    original_rows = len(df)
    output_rows = len(output_df)
    x_added = int(df["X_added"].sum())

    print(f"Input:  {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    print(f"Original genes/rows: {original_rows:,}")
    print(f"Rows after splitting COG categories: {output_rows:,}")
    print(f"New X assignments: {x_added:,}")


if __name__ == "__main__":
    main()
