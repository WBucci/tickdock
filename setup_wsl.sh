#!/bin/bash
PASS="wslP@ss24"

echo "=== [1/4] Install pip + build deps ==="
echo "$PASS" | sudo -S apt-get install -y -qq python3-pip python3-dev cmake 2>&1 | tail -3

echo ""
echo "=== [2/4] Build fpocket from source ==="
cd /tmp
rm -rf fpocket
git clone --depth 1 https://github.com/Discngine/fpocket.git
cd fpocket
# Repo uses a flat Makefile — build in-tree
make -j$(nproc) 2>&1 | grep -E "(Error|fpocket|make\[)" | tail -10
# Binary should be in bin/ or root after successful build
if [ -f "bin/fpocket" ]; then
    echo "$PASS" | sudo -S cp bin/fpocket /usr/local/bin/fpocket
    echo "fpocket installed from bin/"
elif [ -f "fpocket" ]; then
    echo "$PASS" | sudo -S cp fpocket /usr/local/bin/fpocket
    echo "fpocket installed from root"
else
    echo "Binary not found, checking..."
    find . -name "fpocket" -type f 2>/dev/null
    # Try cmake fallback
    mkdir -p cmake_build && cd cmake_build
    cmake -DCMAKE_BUILD_TYPE=Release -S .. -B . 2>&1 | tail -3
    make -j$(nproc) 2>&1 | tail -5
    cd ..
    find . -name "fpocket" -type f 2>/dev/null | head -3
fi
echo "$PASS" | sudo -S chmod 755 /usr/local/bin/fpocket 2>/dev/null || true

echo ""
echo "=== [3/4] Vina already at $(which vina 2>/dev/null || echo NOT_FOUND) ==="
if ! command -v vina &>/dev/null; then
    cd /tmp
    wget -q "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64" -O vina
    chmod +x vina
    echo "$PASS" | sudo -S mv vina /usr/local/bin/vina
fi
echo "vina: $(vina --version 2>&1 | head -1)"

echo ""
echo "=== [4/4] Python packages ==="
python3 -m pip install --quiet --break-system-packages \
    biopython requests pandas rdkit jinja2 tqdm openpyxl
python3 -c "import Bio, requests, pandas, rdkit; print('Python libs: ALL OK')"

echo ""
echo "=== Verification ==="
for tool in fpocket obabel vina; do
    path=$(which $tool 2>/dev/null)
    if [ -n "$path" ]; then echo "  $tool: OK ($path)"
    else echo "  $tool: MISSING"; fi
done

echo ""
echo "=== SETUP COMPLETE ==="
