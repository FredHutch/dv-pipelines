# DynamoDB records for workflows

These are JSON records for the Analysis table. Each folder is a workflow. It's just the JSON record, not the workflow itself. Each file is a workflow for a specific environment (dev, prod, etc).


## Workflow IDs

This is a list of uuids manually created 

**Placeholder Workflows**

```
00000000-0000-0000-0000-000000000000 - 'do nothing' placeholder workflow
11111111-1111-1111-1111-111111111111 - file list placeholder workflow
```

**Cellranger workflows**

```
11111111-2222-1111-1111-111111111111 - Cellranger count
11111111-2222-2222-1111-111111111111 - Cellranger count pubweb
11111111-2222-3333-1111-111111111111 - Cellranger agg pubweb
11111111-2222-3333-1112-111111111111 - Cellranger agg pubweb v2
```

**nf-core workflows**

```
11111111-3333-1111-1111-111111111111 - nf-core scrna-seq
```

**Pubweb workflows**

```
11111111-4444-1111-1111-111111111111 - matrix market to pubweb
11111111-4444-2222-1111-111111111111 - covid-19 atlases to pubweb
```