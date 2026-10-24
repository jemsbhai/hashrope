# Independent controlled-diagnostic review

The completed diagnostic and `controlled_findings.md` pass review. No material
numeric, provenance or interpretation issue remains. This review ran no builds
or timings; it independently read the saved data and recomputed their statistics.

[controlled_review.py](controlled_review.py) verifies the frozen runtime source
copies, archived build inputs, three archived and local executables, saved
context, commands, exit status, sidecars and raw hashes. It independently
recomputes all **248,764 samples, 720 blocks, 36 processes and 90 estimate
vectors**, including the six process ratios, direction counts and log-t
intervals for every metric. All 360 pinned blocks have the required endpoint
processor readings. The machine-readable results are in
[controlled_review.json](controlled_review.json).

The original-executable replication was checked separately: its 24 raw process
outputs, 4,392 samples, 66 paired estimates and retained executable hashes agree
with its summary. Its source archives map exactly to the two controlled cores.
The original and controlled cohorts remain distinct; the identical source
mapping does not imply identical linking, workload wrappers or execution history.

The report accurately retains the independent 1M churn replication estimate
0.9836 [0.9491, 1.0195] and concat estimate 1.0804 [0.9056, 1.2888]. The original
concat point-estimate direction therefore remains unresolved. The four
released/current controlled primary estimates are 0.9622–1.0177 for combined
churn and 1.0014–1.0195 for growing concat. All primary intervals include one;
this is not an equivalence result. The secondary unpinned normal-slot concat
timed-mean ratio 1.0219 [1.0088, 1.0351] and wall ratio 1.0171 [1.0039, 1.0305]
remain visible. Secondary wall intervals also exclude one for an A/A case, which
does not prove that any particular released/current effect is noise.

The report correctly separates slot orientation, affinity modes and A/A controls;
it does not pool them, subtract the separately prepared edit/query timings, or
attribute the historical effect to allocation, processor migration or layout.
The approximately 25 ms calibration target, observed 13.5 ms minimum, 203.8 s
total process wall time and 187 unpinned endpoint changes also match the data.
The no-growth fixture's larger retained arena and different warmup history are
explicit limitations. No causal remedy follows solely from these diagnostics;
the proposed concat optimization requires its own frozen full-workload result.
