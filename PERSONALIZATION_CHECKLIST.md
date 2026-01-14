# Personalization Checklist

Before publishing to GitHub, update these placeholders with your actual information:

## 1. GitHub Username

Replace `ishaq7334` with your actual GitHub ishaq7334 in these files:

### README.md

- Line 110: `git clone https://github.com/ishaq7334/scrnaseq_pipeline.git`
- Line 331: Documentation URL
- Line 384: Citation URL
- Line 405-406: Issues and Discussions links

### mkdocs.yml

- Line 4: `site_url: https://ishaq7334.github.io/scrnaseq_pipeline`
- Line 6: `repo_name: ishaq7334/scrnaseq_pipeline`
- Line 7: `repo_url: https://github.com/ishaq7334/scrnaseq_pipeline`
- Line 87: GitHub profile link
- Line 89: Twitter link (optional)

### docs/index.md

- Line 75: Citation URL
- Line 81-82: Issues and Discussions links
- Line 83: Email

### docs/installation.md

- Line 47: Clone URL

### docs/quickstart.md

- Lines 229-230: Issues and Discussions links

### CONTRIBUTING.md

- Line 59: Clone URL (has `your-ishaq7334` instead of `ishaq7334`)

## 2. Email Address

Replace `ishaaq.raja@gmail.com` with your actual email:

### docs/index.md

- Line 83: `<ishaaq.raja@gmail.com>`

## 3. Quick Update Commands

Run these commands to update everything at once:

```bash
cd /Users/rajaishaqnabikhan/.gemini/antigravity/scratch/scrnaseq_pipeline

# Replace with your GitHub ishaq7334
find . -type f \( -name "*.md" -o -name "*.yml" -o -name "*.yaml" \) \
  -not -path '*/\.git/*' \
  -exec sed -i '' 's/ishaq7334/YOUR_GITHUB_USERNAME/g' {} \;

# Replace with your email
find . -type f -name "*.md" -not -path '*/\.git/*' \
  -exec sed -i '' 's/contact@example\.com/YOUR_EMAIL/g' {} \;

# Commit changes
git add -A
git commit -m "Add personal information"
```

## 4. Example

If your GitHub ishaq7334 is `rajakhan` and email is `raj@example.com`:

```bash
cd /Users/rajaishaqnabikhan/.gemini/antigravity/scratch/scrnaseq_pipeline

find . -type f \( -name "*.md" -o -name "*.yml" -o -name "*.yaml" \) \
  -not -path '*/\.git/*' \
  -exec sed -i '' 's/ishaq7334/rajakhan/g' {} \;

find . -type f -name "*.md" -not -path '*/\.git/*' \
  -exec sed -i '' 's/contact@example\.com/raj@example.com/g' {} \;

git add -A
git commit -m "Add personal information"
```

## 5. Verification

After updating, verify with:

```bash
# Check for remaining placeholders
grep -r "ishaq7334" --include="*.md" --include="*.yml" . | grep -v ".git"
grep -r "ishaaq.raja@gmail.com" --include="*.md" . | grep -v ".git"
```

Should return no results (or only in this checklist file).

## 6. Optional Updates

You may also want to customize:

- **Twitter/Social Links** (mkdocs.yml line 89): Add your Twitter handle or remove
- **Repository Description**: Update in GitHub after creating the repo
- **Topics**: Add relevant topics in GitHub (bioinformatics, single-cell, etc.)

## 7. After Personalization

Once you've updated your information:

1. Test locally: `snakemake --dry-run --cores 1`
2. Commit: `git add -A && git commit -m "Add personal information"`
3. Push to GitHub: `gh repo create scrnaseq_pipeline --public --source=. --push`

---

**Current Status**:

- Author: ✓ RAJ ISHAQ NABI KHAN
- GitHub Username: ⚠️ Needs update (currently: `ishaq7334`)
- Email: ⚠️ Needs update (currently: `ishaaq.raja@gmail.com`)
