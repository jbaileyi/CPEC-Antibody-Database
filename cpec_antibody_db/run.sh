#!/usr/bin/env bash
# One-click launcher for Mac and Linux. Creates a private environment on first run.
cd "$(dirname "$0")"
if [ ! -d .venv ]; then
  echo "First run: installing (about a minute)..."
  python3 -m venv .venv || { echo "Python 3 is required: https://www.python.org/downloads/"; exit 1; }
  .venv/bin/pip install -q -r requirements.txt
fi
.venv/bin/streamlit run app.py
