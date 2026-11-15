#!/bin/bash
set -e
cd "$(dirname "$0")"
echo "EarthTrend Explorer: setting up your local website…"
if ! command -v python3 >/dev/null 2>&1; then
  echo "Python 3 is required. Install it from python.org and run this again."
  read -r -p "Press Enter to close."
  exit 1
fi
if [ ! -d .venv ]; then python3 -m venv .venv; fi
source .venv/bin/activate
python -m pip install -r requirements.txt
(
  for ((attempt=0; attempt<40; attempt++)); do
    if curl -fsS --max-time 1 'http://127.0.0.1:8000/health' >/dev/null 2>&1; then
      open 'http://127.0.0.1:8000'
      break
    fi
    sleep 1
  done
) &
echo "Your website will open at http://127.0.0.1:8000"
echo "To stop the website, press Control+C in this window."
python -m uvicorn backend.app:app --host 127.0.0.1 --port 8000
