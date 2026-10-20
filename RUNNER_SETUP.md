# Self-Hosted GitHub Actions Runner Setup

This guide sets up a self-hosted GitHub Actions runner using Docker Compose with the **myoung34/github-runner** image. rePredRet's model build (`.github/workflows/build-models.yml`) runs on it because a full build can take longer than GitHub-hosted runners allow (6 hours) and needs a persistent cache.

## Runner Scope Options

- **Repository-level**: Only serves one repository (limited use)
- **Organization-level** (what we're setting up): Serves ALL repos in the organization
- **Enterprise-level**: For multiple organizations

**This guide uses organization-level runners**, which is ideal for managing multiple projects. One runner handles jobs from all repos in the organization!

This is equivalent to **GitLab's group-level runners** - just register once at the organization level and it works everywhere.

## Step 0: Prerequisites

- Docker and Docker Compose installed
- Access to `/var/run/docker.sock` on your host
- Local disk space for the runner workspace and the rePredRet cache (the cache holds a clone of rePredRet-models, several GB)
- GitHub account with organization `stanstrup-metabolomics` (already done)

## Step 1: Create a GitHub Personal Access Token (PAT)

This is the easiest authentication method for organization-level runners.

1. Go to: https://github.com/settings/tokens
2. Click **Generate new token** → **Generate new token (classic)**
3. Give it a name: `github-runner`
4. Select scopes:
   - `repo` (full control of private repositories)
   - `admin:org` (required for organization-level runners)
   - `workflow` (for workflow control)
   - `notifications` (optional, for status notifications)
5. Click **Generate token**
6. **Copy the token immediately** (you won't see it again!)

## Step 2: Set Up Environment Variables

Create a `.github_runner.env` file in the directory where you'll run docker-compose:

```bash
# GitHub Personal Access Token (with admin:org scope)
ACCESS_TOKEN=ghp_your_token_here

# Organization name
ORG_NAME=stanstrup-metabolomics

# Runner configuration
RUNNER_NAME_PREFIX=docker-runner
RUNNER_GROUP=default
LABELS=docker,linux,ubuntu
```

Do not set `EPHEMERAL` or `DISABLE_AUTO_UPDATE`, not even to `false`: the image turns
them on when they have any value. An ephemeral runner is deregistered after its first
job and then fails to start with "runner registration has been deleted".

**Note**: The `github_runner.yml` automatically loads `.github_runner.env` into the container - no need to pass `--env-file` in the command!

### Host paths

The host directories mounted into the container are set in a `.env` file next to `github_runner.yml` (Docker Compose reads `.env` for the `${...}` variables in the compose file; `.github_runner.env` only sets variables inside the container). All three are required:

```bash
# Runner configuration and state
RUNNER_DATA_DIR=/path/on/local/disk/github-runner/data
# Runner workspace (/tmp/runner/work in the container; may be cleaned between jobs)
RUNNER_WORK_DIR=/path/on/local/disk/github-runner/work
# Persistent rePredRet cache (/opt/repredret-cache in the container)
REPREDRET_CACHE_DIR=/path/on/local/disk/repredret-cache
```

Use local disks, not network shares: the model build creates ~150k small files and runs git on them, which is very slow over the network.

## Step 3: Start the Runner

```bash
# Navigate to the directory with github_runner.yml, .github_runner.env and .env
cd /path/to/compose

# Start the runner with docker-compose (automatically loads .github_runner.env)
docker-compose --file github_runner.yml up -d

# Check logs (should show "Listening for Jobs")
docker-compose --file github_runner.yml logs -f github-runner
```

**Expected output:**
```
github-runner    | Current runner version: '2.x.x'
github-runner    | Listening for Jobs
```

## Step 4: Verify the Runner is Online and Configure Access

1. Go to https://github.com/organizations/stanstrup-metabolomics/settings/actions/runners
2. You should see your runner with a green dot (online)
3. Click on the runner to open its settings
4. **Important**: Enable "Allow Public Repositories" - this allows the runner to accept jobs from repositories in the organization
5. The runner is now available to ALL repos in the `stanstrup-metabolomics` organization

## Step 5: Workflows Using the Runner

Jobs select the runner by its labels. rePredRet's model build already does:

```yaml
jobs:
  build:
    runs-on: [self-hosted, docker, linux]
```

The same runner handles jobs from **all repos in the organization**. This is equivalent to GitLab's group-level runners!

### What rePredRet's model build needs from the runner

- **`/opt/repredret-cache`** mounted from local disk (`REPREDRET_CACHE_DIR`). The workflow keeps a clone of rePredRet-models in `/opt/repredret-cache/rePredRet-models`, builds models directly in it, and pushes only the changes. The clone is reset to the published state at the start of every run, so a job that died halfway leaves nothing behind.
- **Jobs run directly in the runner container.** The mount is not visible to jobs that use a `container:` image.
- **Git and Git LFS** in the runner image (included in myoung34/github-runner).

## Step 6: Test the Workflow

1. Push a test commit to trigger the workflow, or
2. Go to: https://github.com/stanstrup-metabolomics/rePredRet/actions
3. Click **Run workflow**

You should see the job run on your self-hosted runner.

To test the model build without touching the published models, run **Build and Update Models** manually with `max_systems` set to e.g. 20. Test runs never update rePredRet-models.

## Troubleshooting

### "Listening for Jobs" not appearing in logs
```bash
# Check if token is valid
docker-compose --file github_runner.yml logs github-runner | grep -i error

# Re-create runner with new token
docker-compose --file github_runner.yml down
# Update .github_runner.env with new PAT
docker-compose --file github_runner.yml up -d
```

### Runner not appearing in GitHub UI
- Verify token is valid (hasn't expired)
- Check ORG_NAME matches exactly: `stanstrup-metabolomics`
- Verify token has correct scopes: `repo`, `admin:org`, `workflow`, and `notifications`
- Check logs for authentication errors

### Docker socket permission denied
```bash
# Check socket permissions
ls -la /var/run/docker.sock

# If not readable, add current user to docker group
sudo usermod -aG docker $USER
```

### Check runner logs
```bash
# Real-time logs
docker-compose --file github_runner.yml logs -f github-runner

# Search for specific messages
docker-compose --file github_runner.yml logs | grep -i "registered\|listening\|error"
```

### Runner needs more resources
Edit `github_runner.yml` and increase resource limits:
```yaml
deploy:
  resources:
    limits:
      cpus: '8'
      memory: 16G
    reservations:
      cpus: '4'
      memory: 8G
```

Then restart:
```bash
docker-compose --file github_runner.yml restart
```

## Removing the Runner

1. Stop the Docker container:
```bash
docker-compose --file github_runner.yml down
```

2. (Optional) Delete the host directories set in `.env` (`RUNNER_DATA_DIR`, `RUNNER_WORK_DIR`, `REPREDRET_CACHE_DIR`)

3. Go to https://github.com/organizations/stanstrup-metabolomics/settings/actions/runners
4. Click the three dots next to your runner and select "Remove"

## Running Multiple Runners (Advanced)

To run multiple runners for parallel job execution:

```bash
# Create separate compose files for each runner
cp github_runner.yml github_runner_2.yml

# Create separate env files
cp .github_runner.env .github_runner_2.env

# Edit .github_runner_2.env and update:
# - RUNNER_NAME_PREFIX=docker-runner-2
# - Use different host directories in a separate .env (RUNNER_DATA_DIR, RUNNER_WORK_DIR)
```

Start both:
```bash
docker-compose --file github_runner.yml up -d
docker-compose --file github_runner_2.yml -p runner2 up -d
```

Each runner will be available to the organization and can handle jobs in parallel.

## Notes

- **Docker Image**: `myoung34/github-runner` - actively maintained community image from myoung34
- **Authentication**: Use GitHub Personal Access Token (PAT) with scopes: `repo`, `admin:org`, `workflow`, `notifications`
- **Organization**: `stanstrup-metabolomics` - all runners serve all repos in this org
- **Scope**: `RUNNER_SCOPE: 'org'` makes this runner available to all organization repositories
- **Allow Public Repositories**: Must be enabled in runner settings (https://github.com/organizations/stanstrup-metabolomics/settings/actions/runners) for the runner to accept jobs from organization repositories
- **Persistent mode**: the runner stays registered between jobs as long as `EPHEMERAL` is unset (see Step 2)
- **Storage**: host directories set in `.env` - runner configuration, workspace and the rePredRet cache (local disks)
- **Resource limits**: Adjust based on your system capacity in `github_runner.yml`
- **Docker-in-Docker**: Runner can build and run containers via mounted docker.sock
- **Time limits**: Jobs on self-hosted runners can run for up to 5 days (GitHub-hosted: 6 hours). A job that waits more than 24 hours for an available runner is cancelled, so an offline runner shows up as cancelled runs
- **Inactivity**: GitHub disables scheduled workflows after 60 days without repository activity; re-enable them under Actions
- **Cost**: Free with GitHub Free plan - unlimited runners, only pay for your infrastructure
