#!/usr/bin/env bash
set -e

INPUT="config.h"
OUTPUT="include/tibercad/base/tiber_config.h"

PREFIX="TC_"

mkdir -p "$(dirname "$OUTPUT")"

{
echo "/*"
echo " * This file is auto-generated from config.h"
echo " * Do not edit manually."
echo " */"
echo
echo "#ifndef _TIBER_CONFIG_H_"
echo "#define _TIBER_CONFIG_H_"
echo 
} > "$OUTPUT"

sed -E \
  -e 's/^#undef ([A-Z0-9_]+)/#undef '"$PREFIX"'\1/' \
  -e 's/^#define ([A-Z0-9_]+)/#define '"$PREFIX"'\1/' \
  "$INPUT" >> "$OUTPUT"

echo "" >> "$OUTPUT"
echo "#endif" >> "$OUTPUT"

