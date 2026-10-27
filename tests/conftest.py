"""Run the tests from the repository root, so that the relative data paths of the analysis scripts work."""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
for p in (ROOT, os.path.join(ROOT, "analysis")):
    if p not in sys.path:
        sys.path.insert(0, p)
