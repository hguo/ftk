#!/bin/bash
# analyze_branches.sh - Analyze FTK branches for revitalization planning
#
# Usage: ./scripts/analyze_branches.sh

set -e

echo "=== FTK Branch Analysis ==="
echo "Date: $(date)"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to analyze a branch
analyze_branch() {
    local branch=$1
    local base=${2:-master}

    echo -e "${GREEN}=== Branch: $branch ===${NC}"

    # Check if branch exists
    if ! git rev-parse --verify "$branch" >/dev/null 2>&1; then
        echo -e "${RED}Branch $branch does not exist locally${NC}"
        echo ""
        return
    fi

    # Commits ahead/behind
    local ahead=$(git rev-list --count ${base}..${branch} 2>/dev/null || echo "0")
    local behind=$(git rev-list --count ${branch}..${base} 2>/dev/null || echo "0")
    echo "Commits ahead of $base: $ahead"
    echo "Commits behind $base: $behind"

    # Last commit
    echo "Last commit:"
    git log -1 --oneline "$branch"

    # Changed files
    echo ""
    echo "Changed files (vs $base):"
    git diff --stat "$base...$branch" | tail -1

    # Key file changes
    echo ""
    echo "Key GPU/CUDA files changed:"
    git diff --name-only "$base...$branch" | grep -E '\.(cu|cuh)$' || echo "  None"

    echo ""
    echo "Key header files changed (top 10):"
    git diff --name-only "$base...$branch" | grep -E '\.hh$' | head -10 || echo "  None"

    echo ""
    echo "---"
    echo ""
}

# Analyze priority branches
echo "=== PRIORITY BRANCHES ==="
echo ""

analyze_branch "no_ndarray" "master"
analyze_branch "mpas" "master"
analyze_branch "xgc-velocity" "master"

# Check remote branches
echo "=== REMOTE BRANCHES OF INTEREST ==="
echo ""

for branch in origin/dev-multigpu origin/dev-catch3 origin/dev-union-find-new; do
    if git rev-parse --verify "$branch" >/dev/null 2>&1; then
        analyze_branch "$branch" "master"
    else
        echo -e "${YELLOW}$branch not available locally. Fetch with:${NC}"
        echo "  git fetch origin $(basename $branch)"
        echo ""
    fi
done

# Summary of merge complexity
echo "=== MERGE COMPLEXITY ANALYSIS ==="
echo ""

analyze_merge_complexity() {
    local branch=$1
    local base=${2:-master}

    if ! git rev-parse --verify "$branch" >/dev/null 2>&1; then
        return
    fi

    echo "Merge complexity for $branch -> $base:"

    # Try merge without committing
    git merge-tree $(git merge-base "$base" "$branch") "$base" "$branch" > /tmp/ftk_merge_analysis_$$.txt

    local conflicts=$(grep -c "^changed in both" /tmp/ftk_merge_analysis_$$.txt 2>/dev/null || echo "0")
    local added=$(grep -c "^added in" /tmp/ftk_merge_analysis_$$.txt 2>/dev/null || echo "0")
    local removed=$(grep -c "^removed in" /tmp/ftk_merge_analysis_$$.txt 2>/dev/null || echo "0")

    echo "  Potential conflicts: $conflicts"
    echo "  Files added: $added"
    echo "  Files removed: $removed"

    if [ "$conflicts" -gt 0 ]; then
        echo -e "  ${RED}⚠ Manual conflict resolution needed${NC}"
    else
        echo -e "  ${GREEN}✓ Clean merge expected${NC}"
    fi

    rm -f /tmp/ftk_merge_analysis_$$.txt
    echo ""
}

analyze_merge_complexity "no_ndarray" "master"
analyze_merge_complexity "mpas" "master"
analyze_merge_complexity "xgc-velocity" "master"

# GPU file inventory
echo "=== GPU FILE INVENTORY ==="
echo ""

echo "CUDA files in master:"
find include -name "*.cu" -o -name "*.cuh" 2>/dev/null | wc -l

echo "CUDA files in no_ndarray:"
git ls-tree -r --name-only no_ndarray | grep -E '\.(cu|cuh)$' | wc -l

echo ""
echo "Master CUDA files:"
find include -name "*.cu" -o -name "*.cuh" 2>/dev/null | head -10

echo ""
echo "=== Analysis Complete ==="
