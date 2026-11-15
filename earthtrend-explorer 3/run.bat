@echo off
REM Start Trend Globe on http://localhost:8000   (set DEMO=1 first for synthetic offline data)
cd /d "%~dp0"
if "%DEMO%"=="1" set TREND_GLOBE_DEMO=1
python -m uvicorn backend.app:app --port 8000
