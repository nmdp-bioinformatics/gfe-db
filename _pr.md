Refactor by service using Makefiles

Restructured repository by service:
```bash
.
├── LICENSE
├── Makefile
├── README.md
└── gfe-db
    ├── database             # Database service
    │   ├── Makefile
    │   ├── neo4j
    │   └── template.yaml
    ├── infrastructure       # Network infrastructure including VPC, SSM Parameters and Secrets
    │   ├── Makefile
    │   └── template.yaml
    └── pipeline             # Update pipeline including Batch jobs, StepFunctions, trigger
        ├── Makefile
        ├── config
        ├── functions
        │   ├── Makefile
        │   └── trigger
        ├── jobs
        │   ├── Makefile
        │   ├── build
        │   └── load
        └── template.yaml
```

Tasks accomplished:
- [x] Added local Makefiles to decouple deployment and tear down of individual services and infrastructure
- [x] Root Makefile defines non-sensitive environment variables, secret variables are sourced from `.env` in root
- [x] Root Makefile checks for secret environment variables and cancels deployment if they are not found
- [x] CloudFormation templates use SSM Parameter Store and Secrets Manager so services can access important configuration values faster, easier, and more securely
- [x] EC2 key pair is used if available, created if not during database deployment
- [x] Initial pipeline configuration and repository state is uploaded to S3 after deployment
- [x] Stack deletion occurs in reverse order of deployment to avoid getting stuck from missing dependencies
- [x] Added error handling to Makefiles when deleting resources that don't exist
- [x] DataBucket resource is never deleted, it must be deleted manually (subequent deployments will fail until some logic is added that can handle this condition)