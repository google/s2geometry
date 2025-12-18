# Bazel Central Registry (BCR) Publication

This directory contains configuration files for automated publication of s2geometry to the [Bazel Central Registry (BCR)](https://github.com/bazelbuild/bazel-central-registry).

## Overview

The s2geometry project uses [publish-to-bcr](https://github.com/bazel-contrib/publish-to-bcr) to automatically publish new releases to the Bazel Central Registry. When a new release is created, the GitHub Actions workflow will:

1. Extract release information from the git tag
2. Generate a BCR entry using the template files in this directory
3. Create a pull request to the Bazel Central Registry

## Files in this Directory

### config.yml
Configuration specific to s2geometry:
- `moduleRoots: ["src"]` - Indicates that MODULE.bazel is in the src/ directory

### src/metadata.template.json
Contains the project metadata including:
- Homepage URL
- Maintainer information (Vincent Tsao and Jesse Rosenstock)
- Repository location
- Version tracking (automatically updated)

### src/source.template.json
Defines how to download the source archive:
- `url`: GitHub release archive URL pattern
- `strip_prefix`: Since s2geometry's MODULE.bazel is in the `src/` directory, we strip to `{REPO}-{VERSION}/src`
- `integrity`: Automatically computed by the workflow

### src/presubmit.yml
Defines the BCR presubmit tests that will run when the PR is opened:
- Build targets to verify
- Test targets to run
- Platform matrix (Linux, macOS, macOS ARM64)
- Bazel version matrix (7.x, 8.x, rolling)

## Setup for Maintainers

### Prerequisites

To use the automated BCR publication workflow, a maintainer must:

1. **Create a Personal Access Token (PAT)**
   - Go to GitHub Settings → Developer settings → Personal access tokens → Tokens (classic)
   - Click "Generate new token (classic)"
   - Give it a descriptive name like "BCR Publication for s2geometry"
   - Select the following permissions:
     - ✅ `repo` (Full control of private repositories)
     - ✅ `workflow` (Update GitHub Action workflows)
   - Set an appropriate expiration date (consider using no expiration for automation, or set calendar reminders to regenerate)
   - Click "Generate token" and **copy the token immediately** (you won't be able to see it again)

2. **Add the PAT as a Repository Secret**
   - Go to the s2geometry repository settings
   - Navigate to Secrets and variables → Actions
   - Click "New repository secret"
   - Name: `BCR_PUBLISH_TOKEN`
   - Value: Paste the PAT you created
   - Click "Add secret"

Note: The PAT owner must have a fork of the Bazel Central Registry. The workflow will automatically use the fork associated with the PAT owner's account.

## How to Publish a New Release

Once the setup is complete, publishing to BCR is automatic:

1. **Create a Release on GitHub**
   - Go to the repository's "Releases" page
   - Click "Draft a new release"
   - Create a new tag following the existing pattern (e.g., `v0.13.0`)
   - Fill in release notes
   - Click "Publish release"

2. **Workflow Triggers Automatically**
   - The `publish-to-bcr.yml` workflow will trigger automatically
   - Monitor the workflow run in the Actions tab

3. **Review the BCR Pull Request**
   - Once the workflow completes, a pull request will be opened in the Bazel Central Registry
   - The BCR maintainers will review and merge it
   - You can view the PR by checking the workflow logs for the PR URL

## Manual Triggering

If needed, you can manually trigger the publication workflow:

1. Go to Actions → "Publish to BCR"
2. Click "Run workflow"
3. Enter the tag name (e.g., `v0.13.0`)
4. Click "Run workflow"

## Troubleshooting

### Workflow Fails with "Authentication Error"
- The PAT may have expired or lacks necessary permissions
- Generate a new PAT with `repo` and `workflow` permissions
- Update the `BCR_PUBLISH_TOKEN` secret

### Workflow Fails with "Fork Not Found"
- Ensure the PAT owner has a fork of bazelbuild/bazel-central-registry
- Verify the fork is accessible by the PAT owner

### BCR PR Fails Presubmit Tests
- Check the presubmit.yml configuration in this directory
- Ensure the build and test targets are correct
- The BCR presubmit will show specific error messages

### Module Not Found or Wrong Path
- This usually indicates an issue with `moduleRoots` in config.yml or `strip_prefix` in src/source.template.json
- For s2geometry, MODULE.bazel is in `src/`, so:
  - `config.yml` should have `moduleRoots: ["src"]`
  - `src/source.template.json` should have `strip_prefix: "{REPO}-{VERSION}/src"`
  - Template files are in `.bcr/src/` to match the module structure

## Security Notes

- The PAT should be treated as a sensitive credential
- Only repository administrators should have access to repository secrets
- Consider using a dedicated bot account for BCR publications
- Regularly rotate the PAT (set expiration dates and calendar reminders)
- The PAT only needs write access to your fork of the BCR, not to the main repository

## References

- [publish-to-bcr Documentation](https://github.com/bazel-contrib/publish-to-bcr)
- [Bazel Central Registry](https://github.com/bazelbuild/bazel-central-registry)
- [Bzlmod User Guide](https://bazel.build/docs/bzlmod)
- [s2geometry in BCR](https://github.com/bazelbuild/bazel-central-registry/tree/main/modules/s2geometry)
