#!/bin/sh

# Remove the CluE-* output directories that runs and tests create in the
# repository root, wherever this script is called from.
cd "$(dirname "$0")" || exit 1
rm -rf CluE-*
