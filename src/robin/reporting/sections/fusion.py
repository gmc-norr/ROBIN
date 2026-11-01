"""
fusion.py

This module contains the fusion analysis section of the report.
"""

import os
import logging
import pickle
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.units import inch
from .base import ReportSection
import pandas as pd
import networkx as nx
from reportlab.lib.styles import ParagraphStyle

from robin.classification_config import get_fusion_threshold

logger = logging.getLogger(__name__)

# What the report says for each state of a sample's fusion output.
NOT_RUN_TEXT = "Fusion analysis has not been run for this sample."
UNREADABLE_TEXT = "The fusion results could not be read."
NONE_FOUND_TEXT = "No fusion candidates were identified."
NONE_FOUND_SUMMARY_TEXT = "No fusion networks identified within the target panel."
# When only one candidate type's results are missing and the other's hold candidates.
TARGET_UNREADABLE_SUMMARY_TEXT = "The target-panel fusion results could not be read."
TARGET_NOT_RUN_SUMMARY_TEXT = "Fusion analysis has not been run for the target panel."
GENOME_WIDE_UNREADABLE_SUMMARY_TEXT = "The genome-wide fusion results could not be read."

# The breakpoint validation the report applies on top of the read filters (the same default as
# the GUI's _validate_fusion_breakpoints): reads that support the same breakpoint.
BREAKPOINT_MIN_READS = 4

STATE_NOT_RUN = "not_run"
STATE_UNREADABLE = "unreadable"
STATE_NONE = "none"
STATE_FOUND = "found"


def fusion_methods_sentence():
    """The report's note on how fusion candidates are filtered, built from the config.

    The read filters keep ``mapping_quality > min`` and ``mapping_span > min`` (hence the +1)
    and gene pairs with ``>= read_support`` reads, see fusion_work.
    """
    min_mapq = get_fusion_threshold("mapping_quality") + 1
    min_span = get_fusion_threshold("mapping_span") + 1
    min_reads = get_fusion_threshold("read_support")
    return (
        "Fusion candidates are identified from reads with supplementary alignments and kept "
        f"when both alignments have a mapping quality of at least {min_mapq} and a mapping span "
        f"of at least {min_span} bp, and the gene pair has at least {min_reads} supporting "
        f"reads. Only fusions with {BREAKPOINT_MIN_READS} or more reads at the same validated "
        "breakpoint are reported."
    )


class FusionSection(ReportSection):
    """Section containing the fusion analysis results."""

    def _get_validated_fusion_pairs(self, data):
        """Get validated fusion pairs using the same logic as the GUI.
        
        Uses breakpoint validation with minimum 4 reads support.
        """
        if data is None or data.empty:
            return []
        
        try:
            # Import validation functions from GUI module
            from robin.gui.components.fusion import _cluster_fusion_reads
            
            # Use breakpoint validation (same as GUI)
            clustered_data = _cluster_fusion_reads(
                data, 
                max_distance=10000, 
                use_breakpoint_validation=True
            )
            
            if clustered_data.empty:
                return []
            
            # Extract unique validated gene pairs
            validated_pairs = []
            seen_pairs = set()
            
            for _, row in clustered_data.iterrows():
                fusion_pair_str = row["fusion_pair"]  # e.g., "GENE1-GENE2"
                if fusion_pair_str and fusion_pair_str not in seen_pairs:
                    genes = [g.strip() for g in fusion_pair_str.split("-") if g.strip()]
                    if len(genes) >= 2:
                        # Sort genes for consistency
                        genes_sorted = sorted(genes)
                        pair_key = tuple(genes_sorted)
                        if pair_key not in seen_pairs:
                            validated_pairs.append(pair_key)
                            seen_pairs.add(pair_key)
            
            return validated_pairs
            
        except Exception as e:
            logger.warning(f"Failed to get validated fusion pairs, falling back to simple method: {e}")
            return self._get_gene_pairs_simple(data)

    def _get_gene_pairs_simple(self, data):
        """Simple fallback method for getting gene pairs (without breakpoint validation)."""
        if data is None or data.empty:
            return []

        # For processed data, use read_id grouping
        if 'read_id' in data.columns:
            read_groups = data.groupby("read_id")
        elif 'readID' in data.columns:
            read_groups = data.groupby("readID")
        else:
            return []

        gene_pairs = []
        gene_pair_reads = {}
        
        for _, group in read_groups:
            # Get genes from col4 (Gene column)
            if 'col4' in group.columns:
                genes = group["col4"].unique()
            elif 'Gene' in group.columns:
                genes = group["Gene"].unique()
            else:
                continue
                
            if len(genes) >= 2:
                genes = sorted([str(g).strip() for g in genes if g])
                for i in range(len(genes) - 1):
                    for j in range(i + 1, len(genes)):
                        pair = (genes[i], genes[j])
                        if pair not in gene_pair_reads:
                            gene_pair_reads[pair] = set()
                        # Get read ID
                        if 'read_id' in group.columns:
                            read_id = group["read_id"].iloc[0]
                        elif 'readID' in group.columns:
                            read_id = group["readID"].iloc[0]
                        else:
                            continue
                        gene_pair_reads[pair].add(read_id)

        # Only include pairs with 4 or more supporting reads (matching GUI)
        for pair, reads in gene_pair_reads.items():
            if len(reads) >= 4:
                gene_pairs.append(pair)

        return list(set(gene_pairs))

    def _get_gene_pairs(self, data):
        """Get unique gene fusion pairs from the processed data.
        
        Uses validated fusion pairs with breakpoint validation (same as GUI).
        """
        return self._get_validated_fusion_pairs(data)

    def _get_gene_networks(self, gene_pairs):
        """Create networks of connected genes."""
        G = nx.Graph()
        for pair in gene_pairs:
            G.add_edge(pair[0], pair[1])
        connected_components = list(nx.connected_components(G))
        return [list(component) for component in connected_components]

    def _read_processed_pickle(self, path):
        """Return (candidates, readable) for one processed pickle.

        candidates is the good-pair frame (empty when the pickle holds none). A pickle that
        cannot be read, or holds something else than the processed structure, is not readable.
        """
        try:
            with open(path, "rb") as f:
                processed_data = pickle.load(f)
            annotated_data = processed_data.get("annotated_data", pd.DataFrame())
            goodpairs = processed_data.get("goodpairs", pd.Series(dtype=bool))
            if not isinstance(annotated_data, pd.DataFrame):
                raise TypeError("annotated_data is not a DataFrame")
            # Only keep the good pairs (same as GUI does)
            if not annotated_data.empty and not goodpairs.empty:
                candidates = annotated_data[goodpairs]
            else:
                candidates = annotated_data
            logger.debug(
                f"Loaded fusion candidates from {path} "
                f"({len(candidates)} good pairs from {len(annotated_data)} total candidates)"
            )
            return candidates, True
        except Exception as e:
            logger.warning(f"Could not read fusion pickle {path}: {e}")
            return pd.DataFrame(), False

    def _missing_pickle_state(self, csv_name):
        """State of one candidate type whose processed pickle is missing.

        The type's filtered candidate CSV holds the reads that passed; ``sv_count.txt`` is
        written when the step finishes. CSV rows but no pickle: the results exist but were
        not processed (a folder from before the writer always made one, or this type's build
        failed). No rows: with ``sv_count.txt`` the step ran and nothing passed; without it,
        never run.
        """
        path = os.path.join(self.report.output, csv_name)
        if os.path.exists(path):
            try:
                if not pd.read_csv(path).empty:
                    return STATE_UNREADABLE
            except pd.errors.EmptyDataError:
                pass
            except Exception as e:
                logger.warning(f"Could not read fusion candidates {path}: {e}")
                return STATE_UNREADABLE
        if os.path.exists(os.path.join(self.report.output, "sv_count.txt")):
            return STATE_NONE
        return STATE_NOT_RUN

    def _load_fusion_data(self):
        """Load the processed pickles and decide a state per candidate type and for the sample.

        Per type (``master_state``: target panel; ``all_state``: genome-wide): found when its
        pickle is readable and holds candidates; none when readable and empty; unreadable when
        it can't be read; a missing pickle is decided by that type's own evidence
        (``_missing_pickle_state``). ``state``, the sample's: found if either type is, else
        unreadable if either is, else none if either is, else not run.
        """
        fusion_data = {
            "master_candidates": None,
            "all_candidates": None,
            "master_path": os.path.join(
                self.report.output, "fusion_candidates_master_processed.pkl"
            ),
            "all_path": os.path.join(self.report.output, "fusion_candidates_all_processed.pkl"),
            "master_state": STATE_NOT_RUN,
            "all_state": STATE_NOT_RUN,
            "state": STATE_NOT_RUN,
            "unreadable": [],
        }

        for kind, csv_name in (
            ("master", "fusion_candidates_master.csv"),
            ("all", "fusion_candidates_all.csv"),
        ):
            path = fusion_data[f"{kind}_path"]
            if not os.path.exists(path):
                state = self._missing_pickle_state(csv_name)
            else:
                candidates, ok = self._read_processed_pickle(path)
                if ok:
                    fusion_data[f"{kind}_candidates"] = candidates
                    state = STATE_NONE if candidates.empty else STATE_FOUND
                else:
                    fusion_data["unreadable"].append(path)
                    state = STATE_UNREADABLE
            fusion_data[f"{kind}_state"] = state

        states = {fusion_data["master_state"], fusion_data["all_state"]}
        for state in (STATE_FOUND, STATE_UNREADABLE, STATE_NONE, STATE_NOT_RUN):
            if state in states:
                fusion_data["state"] = state
                break
        return fusion_data

    def _format_fusion_table(self, data, title):
        """Create a formatted table for fusion data using validated fusion pairs."""
        if data is None or data.empty:
            return None

        # Create paragraph styles for cells
        header_style = ParagraphStyle(
            "HeaderStyle",
            parent=self.styles.styles["Normal"],
            fontName="Helvetica-Bold",
            fontSize=8,
            textColor=self.styles.COLORS["primary"],
            leading=10,
            spaceBefore=0,
            spaceAfter=0,
        )

        cell_style = ParagraphStyle(
            "CellStyle",
            parent=self.styles.styles["Normal"],
            fontName="Helvetica",
            fontSize=8,
            textColor=colors.black,
            leading=10,
            spaceBefore=0,
            spaceAfter=0,
        )

        # Prepare table data with updated headers using Paragraph objects
        table_data = [
            [
                Paragraph("Fusion Pair", header_style),
                Paragraph("Chr 1", header_style),
                Paragraph("Chr 2", header_style),
                Paragraph("Gene 1 Breakpoint", header_style),
                Paragraph("Gene 2 Breakpoint", header_style),
                Paragraph("Reads", header_style),
            ]
        ]

        try:
            # Use validated fusion pairs (same as GUI)
            from robin.gui.components.fusion import _cluster_fusion_reads
            
            # Get validated fusion pairs using breakpoint validation
            clustered_data = _cluster_fusion_reads(
                data, 
                max_distance=10000, 
                use_breakpoint_validation=True
            )
            
            if clustered_data.empty:
                return None
            
            # Build fusion details from validated breakpoints
            fusion_details = {}
            for _, row in clustered_data.iterrows():
                fusion_pair_str = row["fusion_pair"]  # e.g., "GENE1-GENE2"
                if fusion_pair_str:
                    genes = [g.strip() for g in fusion_pair_str.split("-") if g.strip()]
                    if len(genes) >= 2:
                        gene_pair = "-".join(sorted(genes))
                        
                        # Get breakpoint information from validated data
                        fusion_details[gene_pair] = {
                            "chrom1": row.get("chr1", "Unknown"),
                            "chrom2": row.get("chr2", "Unknown"),
                            "pos1": row.get("gene1_position", "N/A"),
                            "pos2": row.get("gene2_position", "N/A"),
                            "supporting_reads": int(row.get("reads", 0)),
                        }
            
            # Add rows to table using Paragraph objects for wrappable text
            for gene_pair, details in sorted(fusion_details.items()):
                table_data.append(
                    [
                        Paragraph(gene_pair, cell_style),
                        Paragraph(details["chrom1"], cell_style),
                        Paragraph(details["chrom2"], cell_style),
                        Paragraph(details["pos1"], cell_style),
                        Paragraph(details["pos2"], cell_style),
                        Paragraph(str(details["supporting_reads"]), cell_style),
                    ]
                )

            # Return None if no fusions meet the threshold
            if len(table_data) == 1:  # Only header row
                return None

            # Define column widths (in inches) - adjusted for better proportions
            col_widths = [1.8, 0.5, 0.5, 1.6, 1.6, 0.6]

            # Create and style the table
            table = Table(
                table_data, colWidths=[inch * width for width in col_widths], repeatRows=1
            )
            table.setStyle(
                TableStyle(
                    [
                        # Inherit modern table style
                        *self.MODERN_TABLE_STYLE._cmds,
                        # Preserve specific alignments
                        ("ALIGN", (0, 0), (0, -1), "LEFT"),  # Left-align Fusion Pair
                        ("ALIGN", (1, 0), (2, -1), "CENTER"),  # Center-align Chromosomes
                        ("ALIGN", (3, 0), (4, -1), "LEFT"),  # Left-align positions
                        ("ALIGN", (5, 0), (5, -1), "RIGHT"),  # Right-align Supporting Reads
                    ]
                )
            )

            return table
            
        except Exception as e:
            logger.warning(f"Failed to create fusion table with validation, using fallback: {e}")
            return self._format_fusion_table_fallback(data, title, header_style, cell_style)

    def _format_fusion_table_fallback(self, data, title, header_style, cell_style):
        """Fallback method for formatting fusion table from raw data."""
        # Group by readID to get fusion pairs
        read_groups = data.groupby("readID")

        # Prepare table data with updated headers using Paragraph objects
        table_data = [
            [
                Paragraph("Fusion Pair", header_style),
                Paragraph("Chr 1", header_style),
                Paragraph("Chr 2", header_style),
                Paragraph("Gene 1 Position", header_style),
                Paragraph("Gene 2 Position", header_style),
                Paragraph("Reads", header_style),
            ]
        ]

        # Track processed pairs to avoid duplicates
        processed_pairs = set()
        fusion_details = {}

        # First pass: collect all read IDs for each gene pair
        gene_pair_reads = {}
        for read_id, group in read_groups:
            if len(group) >= 2:
                genes = sorted(
                    group["Gene"].unique()
                )  # Sort to ensure consistent ordering
                if len(genes) >= 2:
                    for i in range(len(genes) - 1):
                        for j in range(i + 1, len(genes)):
                            gene_pair = f"{genes[i]}-{genes[j]}"
                            if gene_pair not in gene_pair_reads:
                                gene_pair_reads[gene_pair] = set()
                            gene_pair_reads[gene_pair].add(read_id)

        # Second pass: process each gene pair with 4 or more supporting reads (matching GUI)
        for gene_pair, reads in gene_pair_reads.items():
            if len(reads) >= 4 and gene_pair not in processed_pairs:
                processed_pairs.add(gene_pair)
                genes = gene_pair.split("-")

                # Get details for each gene in the pair
                gene_info = {}
                for gene in genes:
                    gene_filtered = data[data["Gene"] == gene]
                    if gene_filtered.empty:
                        logger.warning(f"No data found for gene: {gene}")
                        continue
                    gene_data = gene_filtered.iloc[0]
                    gene_info[gene] = {
                        "chrom": gene_data["chromBED"],
                        "position": f"{gene_data['BS']}-{gene_data['BE']}",
                    }

                # Only store fusion details if we have complete information for both genes
                if (
                    len(gene_info) == 2
                    and genes[0] in gene_info
                    and genes[1] in gene_info
                ):
                    fusion_details[gene_pair] = {
                        "chrom1": gene_info[genes[0]]["chrom"],
                        "chrom2": gene_info[genes[1]]["chrom"],
                        "pos1": gene_info[genes[0]]["position"],
                        "pos2": gene_info[genes[1]]["position"],
                        "supporting_reads": len(reads),
                    }
                else:
                    logger.warning(
                        f"Incomplete gene information for fusion pair: {gene_pair}"
                    )

        # Add rows to table using Paragraph objects for wrappable text
        for gene_pair, details in fusion_details.items():
            table_data.append(
                [
                    Paragraph(gene_pair, cell_style),
                    Paragraph(details["chrom1"], cell_style),
                    Paragraph(details["chrom2"], cell_style),
                    Paragraph(details["pos1"], cell_style),
                    Paragraph(details["pos2"], cell_style),
                    Paragraph(str(details["supporting_reads"]), cell_style),
                ]
            )

        # Return None if no fusions meet the threshold
        if len(table_data) == 1:  # Only header row
            return None

        # Define column widths (in inches) - adjusted for better proportions
        col_widths = [1.8, 0.5, 0.5, 1.6, 1.6, 0.6]

        # Create and style the table
        table = Table(
            table_data, colWidths=[inch * width for width in col_widths], repeatRows=1
        )
        table.setStyle(
            TableStyle(
                [
                    # Inherit modern table style
                    *self.MODERN_TABLE_STYLE._cmds,
                    # Preserve specific alignments
                    ("ALIGN", (0, 0), (0, -1), "LEFT"),  # Left-align Fusion Pair
                    ("ALIGN", (1, 0), (2, -1), "CENTER"),  # Center-align Chromosomes
                    ("ALIGN", (3, 0), (4, -1), "LEFT"),  # Left-align positions
                    ("ALIGN", (5, 0), (5, -1), "RIGHT"),  # Right-align Supporting Reads
                ]
            )
        )

        return table

    def add_content(self):
        """Add the fusion analysis content to the report."""
        logger.debug("Starting fusion section content generation")

        # Add section title
        self.elements.append(
            Paragraph("Gene Fusion Analysis", self.styles.styles["Heading1"])
        )
        self.elements.append(Spacer(1, 6))

        # Load fusion data
        fusion_data = self._load_fusion_data()

        # Get unique fusion pairs
        master_pairs = self._get_gene_pairs(fusion_data["master_candidates"])
        all_pairs = self._get_gene_pairs(fusion_data["all_candidates"])

        # Get gene networks
        master_networks = self._get_gene_networks(master_pairs)
        all_networks = self._get_gene_networks(all_pairs)

        # Create summary text: the target-panel sentence follows the target panel's state,
        # the genome-wide sentence the genome-wide state.
        master_state = fusion_data["master_state"]
        summary_text = []
        if master_networks:
            summary_text.append(
                f"Found {len(master_networks)} fusion networks within the target panel:<br/>"
            )
            for network in master_networks:
                summary_text.append(f"• {' - '.join(network)}")
        elif master_state in (STATE_FOUND, STATE_NONE):
            summary_text.append(NONE_FOUND_SUMMARY_TEXT)
        elif fusion_data["state"] in (STATE_FOUND, STATE_NONE):
            # The genome-wide results were read; the target panel's were not.
            summary_text.append(
                TARGET_UNREADABLE_SUMMARY_TEXT
                if master_state == STATE_UNREADABLE
                else TARGET_NOT_RUN_SUMMARY_TEXT
            )
        elif fusion_data["state"] == STATE_UNREADABLE:
            summary_text.append(UNREADABLE_TEXT)
        else:
            summary_text.append(NOT_RUN_TEXT)

        if all_networks and len(all_networks) > len(master_networks):
            summary_text.append(
                f"<br/>Additional {len(all_networks) - len(master_networks)} genome-wide fusion networks identified."
            )
        elif fusion_data["all_state"] == STATE_UNREADABLE and master_state in (
            STATE_FOUND,
            STATE_NONE,
        ):
            # (When the target panel's results are missing too, the sentence above covers both.)
            summary_text.append(f"<br/>{GENOME_WIDE_UNREADABLE_SUMMARY_TEXT}")

        summary_paragraph = "\n".join(summary_text)

        # Add summary to summary section of report
        self.summary_elements.append(
            Paragraph("Gene Fusions", self.styles.styles["Heading3"])
        )
        self.summary_elements.append(
            Paragraph(summary_paragraph, self.styles.styles["Normal"])
        )

        # Add detailed analysis section if there are any candidates
        if master_pairs or all_pairs:
            self.elements.append(
                Paragraph("Detailed Analysis", self.styles.styles["Heading2"])
            )
            self.elements.append(Spacer(1, 6))

            # Add target panel fusions
            if master_pairs:
                self.elements.append(
                    Paragraph("Target Panel Fusions", self.styles.styles["Heading3"])
                )
                table = self._format_fusion_table(
                    fusion_data["master_candidates"], "Target Panel Fusions"
                )
                if table:
                    self.elements.append(table)
                    self.elements.append(Spacer(1, 6))

            # Add genome-wide fusions
            if all_pairs and len(all_pairs) > len(master_pairs):
                self.elements.append(
                    Paragraph("Genome-wide Fusions", self.styles.styles["Heading3"])
                )
                table = self._format_fusion_table(
                    fusion_data["all_candidates"], "Genome-wide Fusions"
                )
                if table:
                    self.elements.append(table)
                    self.elements.append(Spacer(1, 6))

            # Add source information
            source_text = "Data sources:\n"
            source_text += f"- Target panel fusions: {fusion_data['master_path']}\n"
            source_text += f"- Genome-wide fusions: {fusion_data['all_path']}\n"
            source_text += f"\nNote: {fusion_methods_sentence()}"
            if fusion_data["unreadable"]:
                source_text += (
                    "\nThe fusion results in "
                    + ", ".join(fusion_data["unreadable"])
                    + " could not be read."
                )

            self.elements.append(Paragraph(source_text, self.styles.styles["Normal"]))

            # Build export DataFrames
            try:
                import pandas as pd

                # Summaries
                networks_rows = []
                for scope, nets in (
                    ("target", master_networks),
                    ("genome_wide", all_networks),
                ):
                    for net in nets:
                        networks_rows.append(
                            {
                                "Scope": scope,
                                "Network": " - ".join(net),
                                "Genes": ", ".join(net),
                            }
                        )
                if networks_rows:
                    self.export_frames["fusions_networks"] = pd.DataFrame(networks_rows)

                # Detailed tables if available
                def _fusion_rows(df):
                    """Get validated fusion rows for export (same as GUI logic)."""
                    if df is None or df.empty:
                        return []
                    
                    try:
                        # Use validated fusion pairs (same as GUI)
                        from robin.gui.components.fusion import _cluster_fusion_reads
                        
                        # Get validated fusion pairs using breakpoint validation
                        clustered_data = _cluster_fusion_reads(
                            df, 
                            max_distance=10000, 
                            use_breakpoint_validation=True
                        )
                        
                        if clustered_data.empty:
                            return []
                        
                        # Build rows from validated breakpoints
                        rows = []
                        seen_pairs = set()
                        
                        for _, row in clustered_data.iterrows():
                            fusion_pair_str = row["fusion_pair"]  # e.g., "GENE1-GENE2"
                            if fusion_pair_str and fusion_pair_str not in seen_pairs:
                                seen_pairs.add(fusion_pair_str)
                                genes = [g.strip() for g in fusion_pair_str.split("-") if g.strip()]
                                if len(genes) >= 2:
                                    gene_pair = "-".join(sorted(genes))
                                    rows.append({
                                        "FusionPair": gene_pair,
                                        "Chrom1": row.get("chr1", "Unknown"),
                                        "Chrom2": row.get("chr2", "Unknown"),
                                        "Gene1Pos": row.get("gene1_position", "N/A"),
                                        "Gene2Pos": row.get("gene2_position", "N/A"),
                                        "SupportingReads": int(row.get("reads", 0)),
                                    })
                        return rows
                    except Exception as e:
                        logger.warning(f"Failed to get validated fusion rows for export: {e}")
                        # Fallback to simple method with >= 4 reads
                        read_groups = df.groupby("read_id" if "read_id" in df.columns else "readID")
                        gene_pair_reads = {}
                        for read_id, group in read_groups:
                            if 'col4' in group.columns:
                                genes = sorted(group["col4"].unique())
                            elif 'Gene' in group.columns:
                                genes = sorted(group["Gene"].unique())
                            else:
                                continue
                            if len(genes) >= 2:
                                for i in range(len(genes) - 1):
                                    for j in range(i + 1, len(genes)):
                                        gene_pair = f"{genes[i]}-{genes[j]}"
                                        gene_pair_reads.setdefault(gene_pair, set()).add(read_id)
                        rows = []
                        for gene_pair, reads in gene_pair_reads.items():
                            if len(reads) >= 4:  # Minimum 4 reads (matching GUI)
                                # Get gene information
                                if 'col4' in df.columns:
                                    g1_data = df[df['col4'] == gene_pair.split('-')[0]]
                                    g2_data = df[df['col4'] == gene_pair.split('-')[1]]
                                elif 'Gene' in df.columns:
                                    g1_data = df[df["Gene"] == gene_pair.split('-')[0]]
                                    g2_data = df[df["Gene"] == gene_pair.split('-')[1]]
                                else:
                                    continue
                                
                                if not g1_data.empty and not g2_data.empty:
                                    g1row = g1_data.iloc[0]
                                    g2row = g2_data.iloc[0]
                                    
                                    if 'reference_id' in g1row.index:
                                        chrom1 = g1row.get('reference_id', 'Unknown')
                                        chrom2 = g2row.get('reference_id', 'Unknown')
                                        pos1 = f"{g1row.get('reference_start', 'N/A')}-{g1row.get('reference_end', 'N/A')}"
                                        pos2 = f"{g2row.get('reference_start', 'N/A')}-{g2row.get('reference_end', 'N/A')}"
                                    else:
                                        chrom1 = g1row.get('chromBED', 'Unknown')
                                        chrom2 = g2row.get('chromBED', 'Unknown')
                                        pos1 = f"{g1row.get('BS', 'N/A')}-{g1row.get('BE', 'N/A')}"
                                        pos2 = f"{g2row.get('BS', 'N/A')}-{g2row.get('BE', 'N/A')}"
                                    
                                    rows.append({
                                        "FusionPair": gene_pair,
                                        "Chrom1": chrom1,
                                        "Chrom2": chrom2,
                                        "Gene1Pos": pos1,
                                        "Gene2Pos": pos2,
                                        "SupportingReads": len(reads),
                                    })
                        return rows

                master_rows = _fusion_rows(fusion_data["master_candidates"])
                if master_rows:
                    self.export_frames["fusions_target_panel"] = pd.DataFrame(
                        master_rows
                    )

                all_rows = _fusion_rows(fusion_data["all_candidates"])
                if all_rows:
                    self.export_frames["fusions_genome_wide"] = pd.DataFrame(all_rows)
            except Exception:
                pass
        else:
            body_text = {
                STATE_NOT_RUN: NOT_RUN_TEXT,
                STATE_UNREADABLE: UNREADABLE_TEXT,
            }.get(fusion_data["state"], NONE_FOUND_TEXT)
            self.elements.append(
                Paragraph(
                    body_text,
                    self.styles.styles["Normal"],
                )
            )
