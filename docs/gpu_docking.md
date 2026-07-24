# GPU Docking — Vina-GPU 2.1 on AMD RX 9070 XT (RDNA4 / gfx1201)

GPU-accelerated docking for TickDock using **Vina-GPU 2.1** (OpenCL) on the AMD
RX 9070 XT. Keeps AutoDock Vina scoring (consistent with the existing CPU
`top_hits.json` / `pruned_nonhits.jsonl`). **~5× faster than the 16-core CPU**
(measured ~3 ligands/s vs ~0.56/s effective single-target... see Benchmarks).

> Runs **Windows-native** (WSL2 cannot see the AMD GPU). `gpu_screen.py` runs in
> WSL and launches the Windows exe via interop — the exe gets full GPU access.

## Locations
- **Build dir (outside repo):** `%USERPROFILE%\gpu_docking\Vina-GPU-2.1\AutoDock-Vina-GPU-2.1` (set `VINA_GPU_DIR` to override)
- **Binary:** `Vina-GPU-AMD.exe` (+ 6 MinGW DLLs copied alongside → standalone)
- **Cached kernels:** `Kernel1_Opt.bin`, `Kernel2_Opt.bin` (compiled once for gfx1201)
- **Runner (in repo):** `scripts/gpu_screen.py`
- **Build scripts + patch record:** `tools/gpu_build/`

## Toolchain (Windows)
- AMD Adrenalin driver → OpenCL 2.0 runtime (`C:\Windows\System32\OpenCL.dll`), gfx1201
- MSYS2 + MinGW-w64: `pacman -S mingw-w64-x86_64-{gcc,boost,opencl-headers,opencl-icd}`
- (Visual Studio Build Tools NOT needed — the committed Makefile is gcc-based)

## Source patches (vs upstream DeltaGroupNJUPT/Vina-GPU-2.1)
All in `AutoDock-Vina-GPU-2.1/`. See `tools/gpu_build/patches.md` for exact diffs.
1. **`main/main.cpp`** — `boost/filesystem/convenience.hpp` → `boost/filesystem.hpp` (removed in modern Boost)
2. **`lib/main_procedure_cl.cpp`** — add `#include <windows.h>` (for `Sleep`); convert kernel setup to a **runtime cache** (load `Kernel*_Opt.bin` if present, else compile-from-source + save)
3. **`OpenCL/src/wrapcl.cpp`** —
   - build options: `-Werror` → `-w` (AMD's compiler errors on trailing null chars NVIDIA ignored)
   - `SaveProgramToBinary`/`SetupBuildProgramWithBinary`: `fopen "w"/"r"` → **`"wb"/"rb"`** (text mode corrupted the CL binary on Windows → `CL_INVALID_BINARY`)
4. Build flags: `-DAMD_PLATFORM` (not NVIDIA), `-DOPENCL_2_0`, `-DCL_TARGET_OPENCL_VERSION=200`

## Build
```bash
# from MSYS2 (paths are msys-style):
bash $HOME/gpu_docking/restore_and_pack.sh   # builds + copies DLLs
```
First docking run compiles the kernel for gfx1201 (~13s) and writes the `*_Opt.bin`
cache; subsequent runs load it instantly.

## Usage (from WSL)
```bash
python3 scripts/gpu_screen.py --targets B7QK46 --dry-run   # show real gap (valid ligands)
python3 scripts/gpu_screen.py --targets B7P5E9 B7PY20      # screen specific targets
python3 scripts/gpu_screen.py --all                        # all 138 targets
```
Integrates exactly like `fill_target_gaps`: hits → `{target}_results/`, non-hits →
`pruned_nonhits.jsonl`, then rebuilds `top_hits.json`. Keep-awake handled internally.

## Robustness (Vina-GPU is fragile; the runner compensates)
- **Only valid ligands** (`>500B`; skips the empty obabel stubs)
- **Multi-model PDBQTs** (salts, >1 TORSDOF): extract the **largest fragment** (drug, not counterion) — Vina-GPU's parser rejects multi-model
- **Non-standard atom types** (metals, e.g. `Al`): **skipped** — they crash the whole batch (`tree.h:235`), not just the ligand
- **Chunked (512) + bisect-on-crash**: a ligand that still crashes Vina-GPU is isolated via binary-halving and skipped (logged), so one bad ligand never kills the run
- **Windows path quirk:** configs use `C:/...` paths (the exe is native Windows; MSYS `/c/...` fails)

## Benchmarks (B7P5E9, 200 ligands, gfx1201)
Thread sweep (lower = faster, scores barely move):

| thread | throughput | mean score |
|--------|-----------|-----------|
| 1000 | 2.67/s | −9.03 |
| 2000 | ~2.6/s | ~−9.06 |  (default — best balance)
| 8000 | 2.13/s | −9.11 |
| 16000 | 1.61/s | −9.13 |

2-way process concurrency: +16% only (one process ~saturates the GPU). Net ceiling ~3/s.

## Throughput estimates (thread=2000)
| Job | CPU 16-core | GPU |
|-----|------------|-----|
| Finish 5,148 grid (≈342k pairs) | 7–10 days | ~1.3 days |
| 10k × 138 screen (≈1.38M pairs) | ~28 days | ~5 days |

## Methodology note (for the paper)
Vina-GPU uses heuristic `search_depth`, not CPU `exhaustiveness`. Same scoring
function, comparable affinities, but not identical poses. For a clean dataset,
re-dock consistently on GPU rather than mixing GPU + CPU scores.
