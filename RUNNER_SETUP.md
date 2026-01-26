# Self-Hosted GitHub Actions Runner Setup

This guide sets up a self-hosted GitHub Actions runner using Docker Compose to bypass GitHub's job time limits.

## Runner Scope Options

- **Repository-level**: Only serves one repository (limited use)
- **Organization-level** (what we're setting up): Serves ALL repos in the organization
- **Enterprise-level**: For multiple organizations

**This guide uses organization-level runners**, which is ideal for managing multiple projects. One runner handles jobs from all repos in the organization!

This is equivalent to **GitLab's group-level runners** - just register once at the organization level and it works everywhere.

## Step 0: Create the Organization

1. Go to: https://github.com/organizations/new
2. Organization name: `stanstrup-metabolomics`
3. Billing email: (your email)
4. Organization account: `stanstrup-metabolomics`
5. This organization will be: `My personal organization`
6. Click **Create organization**
7. Choose the **Free** plan

## Step 1: Get the Runner Token (Organization Level)

1. Go to your **organization settings**: https://github.com/organizations/stanstrup-metabolomics/settings/actions/runners
   - **NOT** the repository or personal account settings
   - This makes the runner available to ALL repos in the organization

2. Click **New self-hosted runner**
3. Note the token shown in the `Configure` section (valid for 1 hour)

## Step 2: Set Up Environment Variables

Create a `.github_runner.env` file in the rePredRet directory:

```bash
# GitHub repository configuration
REPO_URL=https://github.com/stanstrup/rePredRet

# Runner token (get from Step 1 - only valid for 1 hour)
RUNNER_TOKEN=<paste_token_here>

# Optional: Give the runner a name
RUNNER_NAME=docker-runner-1

# Optional: Labels for the runner (comma-separated)
RUNNER_LABELS=docker,linux,long-jobs

# Optional: Make runner ephemeral (auto-removes after job)
EPHEMERAL=false
```

**Note**: The `github_runner.yml` automatically loads `.github_runner.env` - no need to pass `--env-file` in the command!

## Step 3: Start the Runner

```bash
# Navigate to rePredRet directory
cd /path/to/rePredRet

# Start the runner with docker-compose (automatically loads .github_runner.env)
docker-compose -f github_runner.yml up -d

# Check logs
docker-compose -f github_runner.yml logs -f github-runner
```

## Step 4: Verify the Runner is Online

1. Go to https://github.com/organizations/stanstrup-metabolomics/settings/actions/runners
2. You should see your runner with a green dot (online)
3. This runner will now be available to ALL repos in the `stanstrup-metabolomics` organization

## Step 5 (Optional): Transfer Repos to Organization

To use the runner with your repositories, you can either:

**Option A: Transfer repos to the organization** (recommended)
1. Go to rePredRet repo: https://github.com/stanstrup/rePredRet/settings
2. Scroll down to "Transfer ownership"
3. Select `stanstrup-metabolomics` organization
4. Confirm transfer

Repeat for any other repos you want to use with this runner.

**Option B: Keep repos personal, grant access**
- You can grant the organization access to personal repos if you prefer

## Step 6: Update Workflows to Use Your Runner

In **any repository in the `stanstrup-metabolomics` organization**, edit `.github/workflows/build-models.yml`:

```yaml
jobs:
  build:
    # Change from 'ubuntu-latest' to use your self-hosted runner
    runs-on: [self-hosted, docker, long-jobs]
```

Or use the runner name directly:

```yaml
    runs-on: docker-runner-1
```

The same runner handles jobs from **all repos in the organization**. This is equivalent to GitLab's group-level runners!

## Step 7: Test the Workflow

Push a test commit to trigger the workflow, or manually trigger it from:
https://github.com/stanstrup-metabolomics/rePredRet/actions

## Troubleshooting

### Runner stuck in offline/stale state
```bash
# Stop and remove the runner
docker-compose -f github_runner.yml down

# Remove the container completely
docker rm github-actions-runner

# Get a new token and start again
```

### Check runner logs
```bash
docker-compose -f github_runner.yml logs -f github-runner
```

### Runner needs more resources
Edit `github_runner.yml` and increase the resource limits:
```yaml
deploy:
  resources:
    limits:
      cpus: '8'      # Change from 4 to 8
      memory: 16G    # Change from 8G to 16G
```

### Docker in Docker issues
Make sure the Docker socket is properly mounted:
```bash
ls -la /var/run/docker.sock
```

Should output something like:
```
srw-rw---- 1 root docker /var/run/docker.sock
```

## Removing the Runner

1. Stop the Docker container:
```bash
docker-compose -f github_runner.yml down
```

2. Go to https://github.com/organizations/stanstrup-metabolomics/settings/actions/runners
3. Click the three dots next to your runner and select "Remove"

Note: This removes the runner from ALL repos in the `stanstrup-metabolomics` organization

## Running Multiple Runners (Advanced)

If you want to run multiple jobs in parallel, you can register multiple runners at the organization level:

```bash
# Start runner 1 (uses .github_runner.env by default)
docker-compose -f github_runner.yml up -d

# Start runner 2 with a different config file
# First create .github_runner2.env with a different token and name
# Then override the env_file for this instance:
docker-compose -f github_runner.yml --env-file .github_runner2.env -p runner2 up -d
```

For each runner:
1. Get a NEW token from https://github.com/organizations/stanstrup-metabolomics/settings/actions/runners
2. Create a new `.github_runner2.env`, `.github_runner3.env`, etc.
3. Give each a different `RUNNER_NAME` (e.g., docker-runner-1, docker-runner-2)
4. Start with different project names: default for first, `-p runner2`, `-p runner3`, etc.

Then all runners are available at the organization level and can handle jobs from all repos in parallel!

## Notes

- **Organization name**: `stanstrup-metabolomics` - use this in all organization URLs and when granting access
- **Token expiration**: Runner tokens expire after 1 hour. If your runner can't register, get a new token.
- **Persistent runner**: Set `EPHEMERAL=false` to keep the runner available for multiple jobs (recommended)
- **Ephemeral runner**: Set `EPHEMERAL=true` to auto-remove the runner after each job (more secure, needs re-registration)
- **Resource limits**: Adjust `cpus` and `memory` in `github_runner.yml` based on your machine's capacity
- **Docker socket**: The runner needs access to Docker to build/run containers in workflow steps
- **Free plan**: The organization is on GitHub's free plan, which includes unlimited runners
