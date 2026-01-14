# GitHub Repository Setup Guide

This guide will help you push the scRNA-seq pipeline to GitHub and set up the repository properly.

## Step 1: Create GitHub Repository

### Option A: Using GitHub CLI (Recommended)

```bash
# Install GitHub CLI if not already installed
# macOS: brew install gh

# Authenticate
gh auth login

# Create repository
cd /path/to/scrnaseq_pipeline
gh repo create scrnaseq_pipeline --public --source=. --remote=origin --push

# This will:
# - Create a new public repository on GitHub
# - Add it as remote 'origin'
# - Push your local commits
```

### Option B: Using GitHub Website

1. Go to <https://github.com/new>
2. Repository name: `scrnaseq_pipeline`
3. Description: "Comprehensive single-cell RNA-seq analysis pipeline with Snakemake"
4. Choose Public or Private
5. Do NOT initialize with README (we already have one)
6. Click "Create repository"

Then push your local repository:

```bash
cd /path/to/scrnaseq_pipeline

# Add remote
git remote add origin https://github.com/YOUR_USERNAME/scrnaseq_pipeline.git

# Push to GitHub
git branch -M main
git push -u origin main
```

## Step 2: Configure Repository Settings

### Enable GitHub Pages (for Documentation)

1. Go to repository Settings → Pages
2. Source: Deploy from a branch
3. Branch: `gh-pages` (will be created by CI/CD)
4. Click Save

The documentation will be available at:
`https://YOUR_USERNAME.github.io/scrnaseq_pipeline`

### Add Repository Topics

Add these topics to help others discover your repository:

- `bioinformatics`
- `single-cell`
- `rna-seq`
- `snakemake`
- `seurat`
- `scanpy`
- `genomics`
- `pipeline`
- `workflow`

### Set Repository Description

```
Comprehensive single-cell RNA-seq analysis pipeline with Snakemake, supporting Seurat and Scanpy frameworks
```

### Add Repository URL

If you have a documentation site, add it as the website URL.

## Step 3: Create Initial Release

```bash
# Tag the initial release
git tag -a v1.0.0 -m "Initial release: Complete scRNA-seq pipeline"

# Push tags
git push origin --tags
```

Then on GitHub:

1. Go to Releases → Create a new release
2. Choose tag: v1.0.0
3. Release title: "v1.0.0 - Initial Release"
4. Description: Copy from CHANGELOG.md
5. Publish release

## Step 4: Set Up Branch Protection (Optional)

For collaborative development:

1. Go to Settings → Branches
2. Add rule for `main` branch:
   - Require pull request reviews
   - Require status checks to pass (CI)
   - Require branches to be up to date

## Step 5: Enable GitHub Actions

The CI/CD workflow (`.github/workflows/ci.yml`) will automatically:

- Run linting checks
- Execute tests
- Build documentation
- Deploy to GitHub Pages

First push will trigger the workflow.

## Step 6: Add Collaborators (Optional)

1. Go to Settings → Collaborators
2. Add team members with appropriate permissions

## Step 7: Create Issue Templates

Create `.github/ISSUE_TEMPLATE/bug_report.md`:

```markdown
---
name: Bug report
about: Create a report to help us improve
title: '[BUG] '
labels: bug
assignees: ''
---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Run command '...'
2. With config '...'
3. See error

**Expected behavior**
What you expected to happen.

**Environment:**
 - OS: [e.g. Ubuntu 20.04]
 - Conda version:
 - Snakemake version:

**Additional context**
Add any other context about the problem here.
```

## Step 8: Update README with Correct URLs

After creating the repository, update these placeholders in README.md:

```bash
# Replace YOUR_USERNAME with your actual GitHub username
sed -i '' 's/username/YOUR_ACTUAL_USERNAME/g' README.md
sed -i '' 's/contact@example.com/YOUR_EMAIL/g' README.md

# Commit changes
git add README.md
git commit -m "Update URLs with actual GitHub username"
git push
```

## Step 9: Verify Everything Works

1. Check that repository is visible on GitHub
2. Verify CI/CD workflow runs successfully
3. Confirm documentation builds and deploys
4. Test cloning the repository
5. Ensure all links work

## Step 10: Promote Your Repository

### Add to Lists

- Awesome Bioinformatics lists
- Snakemake workflows catalog
- Seurat/Scanpy community resources

### Share

- Twitter/X with relevant hashtags
- Bioinformatics forums
- Lab/department mailing lists

### Get Feedback

- Ask colleagues to review
- Request feedback on GitHub Discussions
- Present at lab meetings

## Maintenance Checklist

### Regular Tasks

- [ ] Respond to issues within 48 hours
- [ ] Review and merge pull requests
- [ ] Update dependencies quarterly
- [ ] Add new features based on user feedback
- [ ] Keep documentation up to date
- [ ] Tag releases for major updates

### Quality Assurance

- [ ] All tests passing
- [ ] Documentation builds successfully
- [ ] No broken links
- [ ] Code follows style guidelines
- [ ] Changelog updated

## Troubleshooting

### Push Rejected

```bash
# If you get "push rejected" error
git pull origin main --rebase
git push origin main
```

### CI/CD Fails

Check the Actions tab on GitHub for detailed error logs.

### Documentation Not Deploying

1. Ensure GitHub Pages is enabled
2. Check that CI/CD workflow has write permissions
3. Verify `gh-pages` branch exists

## Next Steps

After setting up GitHub:

1. Complete remaining analysis scripts
2. Add example datasets
3. Write comprehensive tutorials
4. Create video walkthroughs
5. Publish preprint/paper
6. Submit to workflow repositories

---

**Repository is now ready for professional use and collaboration!**
