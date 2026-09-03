#!/usr/bin/env bash
set -euo pipefail

if ! command -v cargo >/dev/null 2>&1; then
  echo "Rust is not installed. Install it from https://rustup.rs and rerun this script." >&2
  exit 1
fi

cd "$(dirname "$0")/.."

cargo fmt --all -- --check
cargo check --locked --all-targets --all-features
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-targets --all-features
cargo test --locked --doc
cargo doc --locked --no-deps
cargo build --locked --release
cargo package --locked --allow-dirty

echo "All Virustic2 quality gates passed."
