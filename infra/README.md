# Infrastructure

Terraform config that deploys the demo to AWS:

```
                ┌──────────────┐        ┌──────────────────────┐
   users ─────▶ │  CloudFront  │ ─────▶ │  S3 (frontend bucket)│
                └──────────────┘        └──────────────────────┘

                ┌──────────────┐        ┌──────────────────────┐
   browser ───▶ │ API Gateway  │ ─────▶ │  Lambda (container)  │
                │  (HTTP API)  │        │  FastAPI via Mangum  │
                └──────────────┘        └──────────────────────┘
```

- **Frontend**: Vite build output synced to a private S3 bucket, served through CloudFront via Origin Access Control.
- **Backend**: FastAPI on Lambda as a container image (arm64/Graviton), built from the repo root so `reference/python/mapa.py` and `reference/fixtures/` are included. API Gateway HTTP API proxies requests to it.

## Layout

- `bootstrap/` — one-time config, applied manually with local state. Creates the Terraform remote state backend (S3 + DynamoDB) and the GitHub Actions OIDC deploy role.
- `main.tf`, `variables.tf`, `backend.tf`, `frontend.tf`, `outputs.tf` — application infrastructure, applied by the `deploy` job in `.github/workflows/ci-cd.yml`.

## One-time setup

1. **Apply the bootstrap config** (local state, run once):

   ```bash
   cd infra/bootstrap
   terraform init
   terraform apply
   ```

   Note the three outputs: `state_bucket_name`, `lock_table_name`, `github_actions_role_arn`.

2. **Add GitHub repository secrets** (Settings → Secrets and variables → Actions):

   | Secret | Value |
   |---|---|
   | `AWS_DEPLOY_ROLE_ARN` | `github_actions_role_arn` output |
   | `TF_STATE_BUCKET` | `state_bucket_name` output |
   | `TF_LOCK_TABLE` | `lock_table_name` output |

   Optionally add a repository **variable** `AWS_REGION` if you used a region other than `eu-west-2`.

3. **Push to main** (or trigger the workflow manually via Actions → CI/CD → Run workflow). After the frontend build passes, the deploy job builds the backend image, runs Terraform, builds the frontend against the new API URL, and syncs it to S3.

   The CloudFront URL is printed in the workflow job summary.

## Tearing down

```bash
cd infra
terraform init -backend-config="bucket=<state_bucket_name>" \
  -backend-config="key=mapa-pd-calibration/terraform.tfstate" \
  -backend-config="region=<aws_region>" \
  -backend-config="dynamodb_table=<lock_table_name>"
terraform destroy

cd bootstrap
terraform destroy
```
