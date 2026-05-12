#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SRA Metadata Formatter - Standalone Script
============================================
Format comprehensive SRA metadata into standardized samples summary

This script processes metadata_comprehensive.tsv and generates samples_summary.csv
with standardized fields, SeqType classification, and text cleaning.

Features:
- Field mapping with priority order
- Four-layer SeqType classification strategy
- Text cleaning and standardization
- Synonym mapping
- Quality report generation

Author: Auto-generated from ribo_metadata_workflow
Date: 2025
"""

import os
import re
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from collections import Counter

import pandas as pd


# ============================================================================
# BUILT-IN CONFIGURATION (No external config files needed)
# ============================================================================

# SeqType Classification Keywords
RIBO_KEYWORDS_PRIMARY = [
    'ribo-seq', 'riboseq', 'rpf', 'rfp', 'footprint',
    'ribosome profiling', 'ribosome protected', '80s',
    'translatome', 'translatomic', 'polysome profiling',
    'ribosome footprinting'
]

RIBO_KEYWORDS_CONTEXT = [
    'ribosome', 'polysome', 'protected fragment',
    'monosome', 'disome', 'trisome'
]

RNA_KEYWORDS_PRIMARY = [
    'rna-seq', 'rnaseq', 'transcriptome', 'transcriptomic',
    'total rna', 'mrna', 'mirna', 'lncrna', 'sirna',
    'expression profiling', 'transcriptome sequencing'
]

RNA_KEYWORDS_CONTEXT = [
    'total rna sequencing', 'gene expression',
    'transcriptomic analysis'
]

# Field Mappings (Priority Order)
FIELD_MAPPINGS = {
    'Run': ['Run'],
    'ProjectID': ['BioProject', 'ProjectID', 'project', 'study'],
    'Tissue': [
        'biosample_tissue',
        'biosample_organism_part',
        'biosample_source_name',
        'biosample_body_site',
        'biosample_tissue_type',
        'SampleType',
        'Body_Site'
    ],
    'Treatment': [
        'biosample_treatment',
        'biosample_condition',
        'biosample_growth_condition',
        'biosample_stimulus',
        'biosample_stress',
        'biosample_harvest_info'
    ],
    'Cultivar': [
        'biosample_cultivar',
        'biosample_genotype',
        'biosample_ecotype',
        'biosample_strain',
        'biosample_variety',
        'biosample_line'
    ],
    'Stage': [
        'biosample_developmental_stage',
        'biosample_dev_stage',
        'biosample_age',
        'biosample_growth_stage',
        'biosample_life_stage'
    ],
    'Layout': ['LibraryLayout', 'library_layout', 'Layout'],
    'Species': ['ScientificName', 'Organism', 'organism', 'species'],
    'ReleaseDate': ['ReleaseDate', 'releasedate', 'release_date', 'date'],
    'Replicate': [
        'biosample_replicate',
        'replicate',
        'bio_replicate',
        'Replicate'
    ]
}

# Synonym Mappings
SYNONYMS_TISSUE = {
    'leaf': ['leaves', 'foliage', 'lamina', 'leaf tissue'],
    'root': ['roots', 'root tip', 'root tissue'],
    'seed': ['seeds', 'grain', 'seed tissue'],
    'shoot': ['shoots', 'aerial parts', 'aerial portion'],
    'flower': ['flowers', 'floret', 'inflorescence'],
    'stem': ['stems', 'stem tissue'],
    'whole_plant': ['whole plant', 'entire plant', 'whole seedling']
}

SYNONYMS_TREATMENT = {
    'control': ['ctrl', 'mock', 'untreated', 'wild_type', 'wt'],
    'heat': ['high temperature', 'warm', 'heat stress', 'high_temp'],
    'cold': ['low temperature', 'chill', 'cold stress', 'low_temp', 'chilling'],
    'drought': ['water_deficit', 'dry', 'water_stress', 'dehydration'],
    'salt': ['nacl', 'salinity', 'sodium_chloride'],
    'light': ['light_condition', 'illumination', 'photoperiod']
}

# Layout Standardization
LAYOUT_MAPPING = {
    'SINGLE': 'SE',
    'PAIRED': 'PE',
    'SINGLETON': 'SE'
}

# Text Fields for SeqType Keyword Matching
SEQTYPE_TEXT_FIELDS = [
    'LibraryName',
    'SampleName',
    'Title',
    'title',
    'biosample_title',
    'biosample_sample_title',
    'description',
    'biosample_description',
    'experiment_title'
]


# ============================================================================
# CORE FORMATTING CLASS
# ============================================================================

class SRAMetadataFormatter:
    """Format SRA metadata into standardized samples summary"""

    def __init__(self, verbose: bool = True):
        """
        Initialize formatter

        Args:
            verbose: Print progress information
        """
        self.verbose = verbose

    def determine_seqtype(self, df: pd.DataFrame) -> pd.Series:
        """
        Determine SeqType using four-layer priority strategy

        Args:
            df: Input DataFrame

        Returns:
            Series with SeqType values ('Ribo' or 'RNA')
        """
        seqtypes = []
        total = len(df)

        if self.verbose:
            print(f"[INFO] Determining SeqType for {total} samples...")

        for idx, row in df.iterrows():
            seqtype = 'RNA'  # Default
            confidence = 0.0
            layer_used = 4

            # Layer 1: LibraryStrategy check (Most reliable)
            if 'LibraryStrategy' in df.columns:
                strategy = str(row.get('LibraryStrategy', '')).upper()
                if 'RIBO' in strategy:
                    seqtype = 'Ribo'
                    confidence = 1.0
                    layer_used = 1
                elif any(x in strategy for x in ['RNA', 'TRANSCRIPTOME', 'MRNA', 'MIRNA', 'LNCRNA']):
                    seqtype = 'RNA'
                    confidence = 0.95
                    layer_used = 1

            # Layer 2: BioSample assay_type/seq_type check
            if confidence < 0.95:
                for field in ['biosample_assay_type', 'biosample_seq_type']:
                    if field in df.columns:
                        value = str(row.get(field, '')).lower()
                        if 'ribo' in value or 'footprint' in value or 'translatome' in value:
                            seqtype = 'Ribo'
                            confidence = 0.9
                            layer_used = 2
                            break
                        elif 'rna' in value and 'seq' in value:
                            seqtype = 'RNA'
                            confidence = 0.85
                            layer_used = 2

            # Layer 3: Keyword matching
            if confidence < 0.9:
                # Collect text fields
                text_parts = []
                for field in SEQTYPE_TEXT_FIELDS:
                    if field in df.columns:
                        value = str(row.get(field, ''))
                        if value and value != 'NA':
                            text_parts.append(value.lower())

                combined_text = ' | '.join(text_parts)

                # Calculate keyword scores
                ribo_score = self._calculate_keyword_score(
                    combined_text,
                    RIBO_KEYWORDS_PRIMARY,
                    RIBO_KEYWORDS_CONTEXT
                )
                rna_score = self._calculate_keyword_score(
                    combined_text,
                    RNA_KEYWORDS_PRIMARY,
                    RNA_KEYWORDS_CONTEXT
                )

                # Determine based on scores
                if ribo_score > 0 or rna_score > 0:
                    if ribo_score > rna_score:
                        seqtype = 'Ribo'
                        confidence = (ribo_score - rna_score) / (ribo_score + rna_score + 1)
                        layer_used = 3
                    elif rna_score > ribo_score:
                        seqtype = 'RNA'
                        confidence = (rna_score - ribo_score) / (ribo_score + rna_score + 1)
                        layer_used = 3

            seqtypes.append(seqtype)

        # Print statistics
        if self.verbose:
            seqtype_counts = Counter(seqtypes)
            print(f"[INFO] SeqType distribution:")
            print(f"  - Ribo: {seqtype_counts.get('Ribo', 0)} ({seqtype_counts.get('Ribo', 0)/total*100:.1f}%)")
            print(f"  - RNA: {seqtype_counts.get('RNA', 0)} ({seqtype_counts.get('RNA', 0)/total*100:.1f}%)")

        return pd.Series(seqtypes)

    def _calculate_keyword_score(
        self,
        text: str,
        primary_keywords: List[str],
        context_keywords: List[str]
    ) -> float:
        """
        Calculate keyword matching score

        Args:
            text: Text to search in
            primary_keywords: Primary keywords (weight = 2)
            context_keywords: Context keywords (weight = 1)

        Returns:
            Score (higher is stronger match)
        """
        score = 0.0

        for keyword in primary_keywords:
            if keyword.lower() in text:
                score += 2

        for keyword in context_keywords:
            if keyword.lower() in text:
                score += 1

        return score

    def map_fields(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Map fields from source columns to target fields

        Args:
            df: Input DataFrame

        Returns:
            DataFrame with mapped fields
        """
        if self.verbose:
            print(f"[INFO] Mapping fields for {len(df)} samples...")

        mapped_df = pd.DataFrame()

        for target_field, source_columns in FIELD_MAPPINGS.items():
            # Find first non-empty value
            values = []
            for _, row in df.iterrows():
                value = 'NA'
                for source_col in source_columns:
                    if source_col in df.columns:
                        val = row.get(source_col)
                        if pd.notna(val) and str(val).strip() and str(val) != 'NA':
                            value = str(val).strip()
                            break
                values.append(value)

            mapped_df[target_field] = values

        return mapped_df

    def clean_text(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Clean and standardize text fields

        Args:
            df: Input DataFrame

        Returns:
            DataFrame with cleaned text
        """
        if self.verbose:
            print(f"[INFO] Cleaning and standardizing text...")

        cleaned_df = df.copy()

        for col in cleaned_df.columns:
            if col == 'Species':
                # Keep spaces in species names
                cleaned_df[col] = cleaned_df[col].apply(self._clean_species_name)
            elif col == 'Layout':
                # Standardize layout
                cleaned_df[col] = cleaned_df[col].apply(self._standardize_layout)
            elif col == 'SeqType':
                # Keep SeqType as is
                continue
            else:
                # Standard text cleaning
                cleaned_df[col] = cleaned_df[col].apply(self._clean_field_text)

        return cleaned_df

    def _clean_species_name(self, value: str) -> str:
        """Clean species name (preserve spaces)"""
        if pd.isna(value) or value == 'NA':
            return 'NA'

        # Basic cleanup
        value = str(value).strip()

        # Replace underscores with spaces (scientific names)
        value = value.replace('_', ' ')

        # Remove extra spaces
        value = re.sub(r'\s+', ' ', value)

        return value if value else 'NA'

    def _standardize_layout(self, value: str) -> str:
        """Standardize layout values (SE/PE)"""
        if pd.isna(value) or value == 'NA':
            return 'NA'

        value = str(value).upper().strip()

        # Direct mapping
        if value in LAYOUT_MAPPING:
            return LAYOUT_MAPPING[value]

        # Already standard
        if value in ['SE', 'PE']:
            return value

        # Extract first letter
        if value.startswith('S'):
            return 'SE'
        elif value.startswith('P'):
            return 'PE'

        return 'NA'

    def _clean_field_text(self, value: str) -> str:
        """Clean field text (standardize separators)"""
        if pd.isna(value) or value == 'NA':
            return 'NA'

        value = str(value).strip()

        # Replace spaces and hyphens with underscores
        value = re.sub(r'[\s\-]+', '_', value)

        # Remove special characters at start/end
        value = value.strip('_-')

        # Collapse multiple underscores
        value = re.sub(r'_+', '_', value)

        return value if value else 'NA'

    def apply_synonym_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply synonym mapping for standardization

        Args:
            df: Input DataFrame

        Returns:
            DataFrame with standardized values
        """
        if self.verbose:
            print(f"[INFO] Applying synonym mapping...")

        mapped_df = df.copy()

        # Map tissue synonyms
        if 'Tissue' in mapped_df.columns:
            mapped_df['Tissue'] = mapped_df['Tissue'].apply(
                lambda x: self._map_synonym(str(x).lower(), SYNONYMS_TISSUE)
            )

        # Map treatment synonyms
        if 'Treatment' in mapped_df.columns:
            mapped_df['Treatment'] = mapped_df['Treatment'].apply(
                lambda x: self._map_synonym(str(x).lower(), SYNONYMS_TREATMENT)
            )

        return mapped_df

    def _map_synonym(self, value: str, synonym_dict: Dict[str, List[str]]) -> str:
        """Map value to canonical form using synonym dictionary"""
        if value == 'na' or not value:
            return 'NA'

        value = value.lower().strip()

        # Direct match
        for canonical, synonyms in synonym_dict.items():
            if value == canonical:
                return canonical

        # Synonym match
        for canonical, synonyms in synonym_dict.items():
            if value in [s.lower() for s in synonyms]:
                return canonical

        # Return original (capitalized)
        return value.replace('_', ' ').title().replace(' ', '_')

    def generate_quality_report(self, df: pd.DataFrame) -> str:
        """
        Generate quality report

        Args:
            df: Formatted DataFrame

        Returns:
            Quality report as string
        """
        lines = []
        lines.append("=" * 60)
        lines.append("DATA QUALITY REPORT")
        lines.append("=" * 60)
        lines.append("")

        # Basic statistics
        lines.append(f"Total samples: {len(df)}")
        lines.append("")

        # Field completeness
        lines.append("Field Completeness:")
        for col in df.columns:
            non_na = df[col].apply(lambda x: x != 'NA').sum()
            pct = non_na / len(df) * 100
            lines.append(f"  - {col}: {non_na}/{len(df)} ({pct:.1f}%)")
        lines.append("")

        # SeqType distribution
        if 'SeqType' in df.columns:
            lines.append("SeqType Distribution:")
            seqtype_counts = df['SeqType'].value_counts()
            for seqtype, count in seqtype_counts.items():
                pct = count / len(df) * 100
                lines.append(f"  - {seqtype}: {count} ({pct:.1f}%)")
            lines.append("")

        # Layout distribution
        if 'Layout' in df.columns:
            lines.append("Layout Distribution:")
            layout_counts = df['Layout'].value_counts()
            for layout, count in layout_counts.items():
                pct = count / len(df) * 100
                lines.append(f"  - {layout}: {count} ({pct:.1f}%)")
            lines.append("")

        # Top values for key fields
        for field in ['Tissue', 'Treatment', 'Cultivar', 'Stage']:
            if field in df.columns:
                # Exclude 'NA'
                values = df[df[field] != 'NA'][field]
                if len(values) > 0:
                    lines.append(f"Top {field} values:")
                    value_counts = values.value_counts().head(10)
                    for val, count in value_counts.items():
                        pct = count / len(df) * 100
                        lines.append(f"  - {val}: {count} ({pct:.1f}%)")
                    lines.append("")

        lines.append("=" * 60)

        return "\n".join(lines)

    def format(
        self,
        input_file: Path,
        output_file: Path,
        seqtype_override: Optional[str] = None
    ) -> Dict:
        """
        Main formatting pipeline

        Args:
            input_file: Input TSV file path
            output_file: Output CSV file path
            seqtype_override: Override SeqType for all samples ('Ribo' or 'RNA')

        Returns:
            Dictionary with results and statistics
        """
        # Step 1: Read input
        if self.verbose:
            print("\n" + "=" * 60)
            print("SRA METADATA FORMATTING")
            print("=" * 60)
            print(f"[INFO] Reading input file: {input_file}")

        df = pd.read_csv(input_file, sep='\t', dtype=str)
        print(f"[INFO] Loaded {len(df)} samples")

        # Step 2: Map fields
        print("\n" + "-" * 60)
        mapped_df = self.map_fields(df)

        # Step 3: Determine SeqType
        print("\n" + "-" * 60)
        if seqtype_override:
            print(f"[INFO] Overriding SeqType to: {seqtype_override}")
            mapped_df['SeqType'] = seqtype_override
        else:
            mapped_df['SeqType'] = self.determine_seqtype(df)

        # Step 4: Clean text
        print("\n" + "-" * 60)
        cleaned_df = self.clean_text(mapped_df)

        # Step 5: Apply synonym mapping
        print("\n" + "-" * 60)
        final_df = self.apply_synonym_mapping(cleaned_df)

        # Step 6: Select output columns
        output_columns = [
            'Run', 'SeqType', 'ProjectID', 'Tissue',
            'Treatment', 'Cultivar', 'Stage', 'Layout', 'Species'
        ]

        # Ensure all columns exist
        for col in output_columns:
            if col not in final_df.columns:
                final_df[col] = 'NA'

        output_df = final_df[output_columns]

        # Step 7: Save output
        print("\n" + "-" * 60)
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        output_df.to_csv(output_file, index=False, encoding='utf-8')
        print(f"[INFO] Saved formatted output: {output_file}")

        # Step 8: Generate quality report
        print("\n" + "-" * 60)
        quality_report = self.generate_quality_report(output_df)

        report_file = output_file.parent / (output_file.stem + '_quality_report.txt')
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(quality_report)
        print(f"[INFO] Saved quality report: {report_file}")

        # Print summary
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"Total samples processed: {len(output_df)}")
        print(f"Output fields: {len(output_columns)}")
        print(f"\nOutput file: {output_file}")
        print(f"Quality report: {report_file}")

        return {
            'output_df': output_df,
            'quality_report': quality_report,
            'output_file': str(output_file),
            'report_file': str(report_file)
        }


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Format comprehensive SRA metadata into standardized samples summary",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python format_sra_metadata.py \\
    --input metadata_comprehensive.tsv \\
    --output samples_summary.csv

  # Override SeqType (force all samples as Ribo)
  python format_sra_metadata.py \\
    --input metadata_comprehensive.tsv \\
    --output samples_summary.csv \\
    --seq-type-override Ribo

  # Quiet mode (less verbose)
  python format_sra_metadata.py \\
    --input metadata_comprehensive.tsv \\
    --output samples_summary.csv \\
    --quiet
        """
    )

    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Input TSV file (metadata_comprehensive.tsv from batch_fetch_sra_metadata.py)'
    )

    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output CSV file (samples_summary.csv)'
    )

    parser.add_argument(
        '--seq-type-override',
        type=str,
        choices=['Ribo', 'RNA'],
        default=None,
        help='Override SeqType for all samples (default: auto-detect)'
    )

    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Quiet mode (less verbose output)'
    )

    args = parser.parse_args()

    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {args.input}")
        sys.exit(1)

    # Initialize formatter
    formatter = SRAMetadataFormatter(verbose=not args.quiet)

    # Run formatting
    try:
        result = formatter.format(
            input_file=input_path,
            output_file=Path(args.output),
            seqtype_override=args.seq_type_override
        )
        print("\n[INFO] Done!")
        sys.exit(0)
    except Exception as e:
        print(f"\n[ERROR] Formatting failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
