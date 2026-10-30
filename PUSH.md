# Pushing this repo to GitHub

The repository is already initialised and committed. From this folder:

```bash
gh repo create gateway-cloning --public --source=. --remote=origin --push \
  --description "Offline single-file workbench for Gateway cloning: BP/LR simulation on real plasmid sequences, protocol-based attB primer design, reading-frame checks and GenBank export"
```

Then turn on GitHub Pages so the app is live at
<https://linjianglao.github.io/gateway-cloning/>:

```bash
gh api -X POST repos/linjianglao/gateway-cloning/pages \
  -F 'source[branch]=main' -F 'source[path]=/'
```

Or in the browser: **Settings → Pages → Source: Deploy from a branch → main / (root)**.

Pages takes a minute or two to build the first time.

## Without the gh CLI

Create an empty repo named `gateway-cloning` on github.com (no README, no .gitignore,
no licence — this repo already has them), then:

```bash
git remote add origin https://github.com/linjianglao/gateway-cloning.git
git push -u origin main
```
