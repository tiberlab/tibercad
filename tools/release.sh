#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Usage
# -----------------------------
if [ $# -ne 1 ]; then
  echo "Usage: $0 <version>"
  echo "Example: $0 3.5.0"
  exit 1
fi

VERSION="$1"
TAG="v$VERSION"
VERSION_FILE="VERSION"

# -----------------------------
# Helpers
# -----------------------------
die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_clean_tree() {
  git diff --quiet || die "Working tree is dirty"
  git diff --cached --quiet || die "Index is dirty"
}

# -----------------------------
# Preconditions
# -----------------------------

# Must be in top-level repo
git rev-parse --show-toplevel >/dev/null

# Ensure clean main repo
require_clean_tree

# Ensure submodules are clean and on exact tags
git submodule foreach --quiet '
  git diff --quiet || {
    echo "Submodule $name has uncommitted changes"
    exit 1
  }

  git describe --tags --exact-match >/dev/null || {
    echo "Submodule $name is not on an exact tag"
    exit 1
  }
'

# Ensure tag does not already exist
git rev-parse "$TAG" >/dev/null 2>&1 && die "Tag $TAG already exists"

# -----------------------------
# Write VERSION file
# -----------------------------
echo "$VERSION" > "$VERSION_FILE"

# -----------------------------
# Regenerate autotools
# -----------------------------
autoreconf -fi

# -----------------------------
# Commit release state
# -----------------------------
git add "$VERSION_FILE" configure aclocal.m4 Makefile.in */Makefile.in
git commit -m "Release $TAG"

# -----------------------------
# Create annotated tag
# -----------------------------
git tag -a "$TAG" -m "Release $TAG"

echo
echo "✔ Release $TAG created successfully"
echo
echo "Next steps:"
echo "  git push origin main"
echo "  git push origin $TAG"

