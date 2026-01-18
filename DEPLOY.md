# GitHub Repository Setup and PyPI Publishing Guide

## Step 1: Initialize Git Repository

```bash
cd /Users/fcliquet/Workspace/genetics-viz
git init
git add .
git commit -m "Initial commit: genetics-viz web application"
```

## Step 2: Create GitHub Repository

1. Go to https://github.com/organizations/bourgeron-lab/repositories/new
2. Repository name: `genetics-viz`
3. Description: "A web-based visualization tool for genetics cohort data"
4. Set to **Public** (required for PyPI publishing)
5. **DO NOT** initialize with README, license, or gitignore (we already have them)
6. Click "Create repository"

## Step 3: Push to GitHub

```bash
git remote add origin https://github.com/bourgeron-lab/genetics-viz.git
git branch -M main
git push -u origin main
```

## Step 4: Configure PyPI Trusted Publishing

### 4.1 Go to PyPI and set up Trusted Publishing:

1. Log in to https://pypi.org/
2. Go to your account settings
3. Navigate to "Publishing" section
4. Click "Add a new pending publisher"
5. Fill in:
   - **PyPI Project Name**: `genetics-viz`
   - **Owner**: `bourgeron-lab`
   - **Repository name**: `genetics-viz`
   - **Workflow name**: `publish.yml`
   - **Environment name**: `release`
6. Click "Add"

### 4.2 Create GitHub Environment (for protection):

1. Go to https://github.com/bourgeron-lab/genetics-viz/settings/environments
2. Click "New environment"
3. Name it: `release`
4. (Optional) Add deployment protection rules:
   - Required reviewers
   - Wait timer
5. Save

## Step 5: Create and Push Your First Release

When you're ready to publish version 0.1.0:

```bash
# Make sure all changes are committed
git add .
git commit -m "Release v0.1.0"

# Create and push the tag
git tag v0.1.0
git push origin v0.1.0
```

The GitHub Action will automatically:
1. Build the package
2. Publish it to PyPI

## Step 6: Monitor the Release

1. Go to https://github.com/bourgeron-lab/genetics-viz/actions
2. Watch the "Publish to PyPI" workflow run
3. Once complete, your package will be available at: https://pypi.org/project/genetics-viz/

## Future Releases

For subsequent releases:

1. Update version in `pyproject.toml`
2. Commit changes:
   ```bash
   git add pyproject.toml
   git commit -m "Bump version to 0.2.0"
   ```
3. Create and push tag:
   ```bash
   git tag v0.2.0
   git push origin main
   git push origin v0.2.0
   ```

## Version Numbering

Follow semantic versioning (SemVer):
- **MAJOR** version (1.0.0): Incompatible API changes
- **MINOR** version (0.1.0): Add functionality (backwards-compatible)
- **PATCH** version (0.0.1): Bug fixes (backwards-compatible)

## Troubleshooting

### If the PyPI publish fails:

1. Check the GitHub Actions log for errors
2. Verify the PyPI Trusted Publishing configuration matches exactly
3. Ensure the `release` environment exists in GitHub
4. Make sure the repository is public

### To test locally before publishing:

```bash
uv build
# Check dist/ folder for wheel and tar.gz files
```

## Installing Your Package

Once published, users can install with:

```bash
pip install genetics-viz
# or
uv pip install genetics-viz
```

## Notes

- The MIT LICENSE file has been created
- The .github/workflows/publish.yml uses PyPI's Trusted Publishing (no API tokens needed)
- The pyproject.toml has been updated with bourgeron-lab URLs
- Tags must start with 'v' followed by semantic version (e.g., v0.1.0, v1.2.3)
