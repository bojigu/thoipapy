"""THOIPA: machine-learning prediction of residues driving homotypic transmembrane interactions.

Package __init__ files are empty by house style -- import from the defining module rather than
from the package root. The single exception below is deliberate: run_THOIPA_prediction is the
documented public entry point, and re-exporting it keeps
`import thoipapy; thoipapy.run_THOIPA_prediction(...)` working for existing callers.

Before this, every __init__ imported all of its submodules, so `import thoipapy` transitively
pulled in 70 of the 85 modules -- the whole ML stack, seaborn, statsmodels and a package of
deprecated code -- on every process start, including for a single prediction.
"""
from thoipapy.predict import run_THOIPA_prediction

__all__ = ["run_THOIPA_prediction"]
