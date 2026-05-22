#!/bin/bash
echo "=== Tool check ==="
for t in fpocket obabel vina; do
    path=$(which $t 2>/dev/null)
    if [ -n "$path" ]; then echo "  $t: OK ($path)"
    else echo "  $t: MISSING"; fi
done

echo ""
echo "=== Python packages ==="
python3 -c "import Bio, requests, pandas, rdkit, jinja2, tqdm, openpyxl; print('  ALL OK')"

echo ""
echo "=== Pipeline prereq check ==="
python3 /mnt/c/Users/Owner/Documents/AndroidApps/TTD/run_pipeline.py --check
