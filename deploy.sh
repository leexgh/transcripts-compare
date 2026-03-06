#!/bin/bash
set -e

# Deploy to GitHub Pages
# Usage: ./deploy.sh
# Prerequisites:
#   1. Create a GitHub repo
#   2. Update REPO below with your repo URL
#   3. Enable GitHub Pages in repo settings (source: gh-pages branch)
#   4. Set `base` in vite.config.ts to '/repo-name/' for GitHub Pages subpath

REPO="git@github.com:leexgh/transcript-compare.git"

cd "$(dirname "$0")"
npm run build
mkdir -p dist/data
cp public/data/gene_data.json dist/data/gene_data.json

cd dist
git init
git add -A
git commit -m "Deploy"
git push -f "$REPO" main:gh-pages

echo "Deployed! Site will be available at your GitHub Pages URL shortly."
