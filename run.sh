#!/usr/bin/env sh
# Start Trend Globe on http://localhost:8000   (DEMO=1 ./run.sh for synthetic offline data)
cd "$(dirname "$0")"
if [ "$DEMO" = "1" ]; then export TREND_GLOBE_DEMO=1; fi
python -m uvicorn backend.app:app --port "${PORT:-8000}"
