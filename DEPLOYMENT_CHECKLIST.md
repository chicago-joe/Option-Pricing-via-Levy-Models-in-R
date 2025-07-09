# 🚀 Deployment Checklist

Use this checklist to ensure successful deployment of your MyST site to GitHub Pages.

## Pre-Deployment

- [ ] PDF files are placed in `documentation/` folder:
  - [ ] `Project Research Proposal.pdf`
  - [ ] `Option Pricing in Levy Models - Feng et al - Academic Paper.pdf`
  - [ ] `R-Finance Presentation Slides.pdf`
- [ ] All markdown files are properly formatted
- [ ] Images (if any) are placed in `assets/` folder

## GitHub Setup

- [ ] Create/update repository on GitHub
- [ ] Enable GitHub Pages:
  1. Go to Settings → Pages
  2. Set Source to "GitHub Actions"

## Deployment Steps

1. [ ] Initialize git repository (if needed):
   ```bash
   git init
   ```

2. [ ] Add all files:
   ```bash
   git add .
   ```

3. [ ] Commit changes:
   ```bash
   git commit -m "Deploy MyST site for Options Pricing research"
   ```

4. [ ] Add remote origin (replace with your repo URL):
   ```bash
   git remote add origin https://github.com/chicago-joe/Option-Pricing-via-Levy-Models.git
   ```

5. [ ] Push to GitHub:
   ```bash
   git push -u origin main
   ```

## Post-Deployment

- [ ] Check Actions tab on GitHub for build status
- [ ] Wait for green checkmark indicating successful deployment
- [ ] Visit your site at: `https://chicago-joe.github.io/Option-Pricing-via-Levy-Models/`
- [ ] Test all three pages load correctly
- [ ] Verify PDF download links work

## Troubleshooting

If deployment fails:
- [ ] Check Actions tab for error messages
- [ ] Verify all file paths are correct
- [ ] Ensure PDFs are not corrupted
- [ ] Check that `myst.yml` is valid YAML

## Success! 🎉

Once all items are checked, your MyST site should be live and accessible!
