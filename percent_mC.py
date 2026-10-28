#!/usr/bin/env python3
"""Kept so older commands still work. Same as: python3 methyl_aligner.py plot ..."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import methyl_aligner

if __name__ == '__main__':
    methyl_aligner.main(['plot'] + sys.argv[1:])
