#!/bin/sh
# Extensive TMD phonon study driver.
# 6 materials in parallel (20 threads each = 120 of 128 cores), then a
# supercell-convergence check on MoS2 and WSe2 (reps 3 -> 5).
cd "$HOME/QATK_RUNS/tmd" || exit 1
ATK=/home/quantumatk/quantumatk/Y-2026.03/bin/atkpython

echo "=== BATCH 1: TMD family, reps=3, 6 jobs x 20 threads ==="
for M in MoS2 WS2 MoSe2 WSe2 MoTe2 WTe2; do
    $ATK --threads 20 tmd_extensive.py "$M" 3 > "ext_${M}_r3.log" 2>&1 &
done
wait
echo "=== BATCH 1 DONE ==="

echo "=== BATCH 2: supercell convergence, reps=5, 2 jobs x 60 threads ==="
for M in MoS2 WSe2; do
    $ATK --threads 60 tmd_extensive.py "$M" 5 > "ext_${M}_r5.log" 2>&1 &
done
wait
echo "=== BATCH 2 DONE ==="

echo "ALL_EXTENSIVE_DONE"
