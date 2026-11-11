@echo off
cd /d "%~dp0"
if not exist .venv (
  echo First run: installing, about a minute...
  python -m venv .venv || (echo Python 3 is required: https://www.python.org/downloads/ & pause & exit /b 1)
  .venv\Scripts\pip install -q -r requirements.txt
)
.venv\Scripts\streamlit run app.py
pause
