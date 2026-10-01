"""limen: where a molecular distinction becomes visible to sequencing."""
from .identify import (Gene, GeneResult, Transcript, build_matrix,
                       identifiability, incidence_patterns, private_fraction,
                       read_signature, signature_owners, sweep)
from .annotation import load_annotation, detect_format

__version__ = "0.1.0"
__all__ = ["Gene", "Transcript", "GeneResult", "identifiability", "sweep",
           "build_matrix", "incidence_patterns", "private_fraction",
           "read_signature", "signature_owners",
           "load_annotation", "detect_format", "__version__"]
