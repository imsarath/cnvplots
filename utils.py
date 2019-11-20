import pandas as pd


def _get_hover_text(df, genename=None, position=None, log2=None):
    """
    CNV plot - hover text
    """
    hover_text = ''

    if genename is not None and genename in df.columns:
        hover_text = hover_text + '<br>GENE: ' \
                                + df[genename].astype(str)

    if position is not None and position in df.columns:
        hover_text = hover_text + '<br>POSITION: ' \
                                + df[position].astype(str)

    if log2 is not None and log2 in df.columns:
        hover_text = hover_text + '<br>LOG2: ' \
                                + df[log2].astype(str)

    return hover_text


def by_chromosome(segments):

    for chrom, subtable in segments.groupby("chromosome", sort=False):
        yield chrom, pd.DataFrame(subtable)
