variable "aws_region" {
  description = "AWS region for the state bucket, lock table and IAM resources."
  type        = string
  default     = "eu-west-2"
}

variable "project_name" {
  description = "Short name used as a prefix for all resources."
  type        = string
  default     = "mapa-pd-calibration"
}

variable "github_repository" {
  description = "GitHub repository allowed to assume the deploy role, as \"owner/repo\"."
  type        = string
  default     = "dillonsnyman1/mapa-pd-calibration"
}
