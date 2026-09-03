# Contributing to Virustic2

Focused, test-backed changes are welcome. Before opening a pull request, run:

```bash
./scripts/verify.sh
```

Algorithm changes must state the biological assumption they introduce and include a fixture where
that assumption should *not* apply. Parser changes need malformed-input coverage. Performance claims
must include the dataset digest, exact command, toolchain, machine, wall time, and peak resident memory.

Do not silently weaken paired-read validation, strand symmetry, complete edge traversal, deterministic
output, or atomic file behavior. JSON schema changes require a changelog entry.
