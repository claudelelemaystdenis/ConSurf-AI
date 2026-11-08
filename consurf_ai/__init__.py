"""ConSurf-AI: per-residue evolutionary conservation grades without a multiple sequence alignment."""
__version__ = "1.0.0"


def __getattr__(name):            # lazy import so `import consurf_ai` does not load torch/transformers
    if name == "ConSurfAI":
        from .predictor import ConSurfAI
        return ConSurfAI
    raise AttributeError(name)
