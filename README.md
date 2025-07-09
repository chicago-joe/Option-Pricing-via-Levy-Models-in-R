# Options Pricing in Lévy Models - MyST Site

This repository contains a MyST (Markedly Structured Text) multi-page documentation site for the Options Pricing in Lévy Models research project.

## 📁 Project Structure

```
option-pricing-myst-site/
├── .github/
│   └── workflows/
│       └── deploy-myst.yml      # GitHub Actions workflow for deployment
├── _static/
│   └── custom.css              # Custom styling for the site
├── assets/                     # Images and other assets
├── documentation/              # PDF documents
│   ├── Project Research Proposal.pdf
│   ├── Option Pricing in Levy Models - Feng et al - Academic Paper.pdf
│   └── R-Finance Presentation Slides.pdf
├── myst.yml                    # MyST configuration file
├── index.md                    # Main page: Research Proposal
├── academic-paper.md           # Academic Paper page
├── presentation-slides.md      # R/Finance Presentation page
└── README.md                   # This file
```

## 🚀 Deployment Instructions

### Prerequisites

1. **GitHub Repository**: Ensure your code is pushed to GitHub
2. **GitHub Pages**: Enable GitHub Pages in your repository settings
3. **PDFs**: Place your PDF files in the `documentation/` folder

### Step 1: Enable GitHub Pages

1. Go to your repository on GitHub
2. Navigate to **Settings** → **Pages**
3. Under **Source**, select **GitHub Actions**

### Step 2: Add PDF Files

Place your three PDF files in the `documentation/` directory:
- `Project Research Proposal.pdf`
- `Option Pricing in Levy Models - Feng et al - Academic Paper.pdf`
- `R-Finance Presentation Slides.pdf`

### Step 3: Push to GitHub

```bash
# Initialize git (if not already done)
git init

# Add all files
git add .

# Commit
git commit -m "Initial MyST site setup for Options Pricing research"

# Add your GitHub repository as origin
git remote add origin https://github.com/chicago-joe/Option-Pricing-via-Levy-Models.git

# Push to main branch
git push -u origin main
```

### Step 4: Monitor Deployment

1. Go to the **Actions** tab in your GitHub repository
2. You should see the "Deploy MyST Site" workflow running
3. Once complete (green checkmark), your site will be available at:
   ```
   https://chicago-joe.github.io/Option-Pricing-via-Levy-Models/
   ```

## 📝 Content Management

### Adding/Editing Pages

1. Edit the markdown files (`index.md`, `academic-paper.md`, `presentation-slides.md`)
2. Update `myst.yml` if adding new pages
3. Push changes to trigger automatic deployment

### Adding Images

1. Place images in the `assets/` directory
2. Reference in markdown using:
   ```markdown
   ```{figure} assets/your-image.png
   :width: 600px
   :align: center
   :alt: Image description
   
   Caption for the image
   ```

### Updating PDFs

Simply replace the PDF files in the `documentation/` directory and push the changes.

## 🎨 Customization

### Styling

Edit `_static/custom.css` to customize the appearance of your site.

### Site Configuration

Modify `myst.yml` to update:
- Site title and metadata
- Navigation structure
- Author information
- Export settings

## 🛠️ Local Development

To test the site locally:

```bash
# Install MyST
pip install mystmd

# Build the site
myst build --html

# Serve locally
myst start
```

The site will be available at `http://localhost:3000`

## 📊 Features

- **Multi-page structure**: Three interconnected pages
- **PDF integration**: Direct links to research documents
- **Responsive design**: Works on all devices
- **Code highlighting**: R code examples with syntax highlighting
- **Mathematical equations**: LaTeX support for formulas
- **Interactive elements**: Mermaid diagrams and tables
- **Dark mode**: Automatic dark mode support

## 🐛 Troubleshooting

### Build Failures

1. Check the Actions tab for error messages
2. Ensure all markdown files are valid
3. Verify PDF files exist in the correct location

### Page Not Found (404)

1. Wait a few minutes after deployment
2. Check that GitHub Pages is enabled
3. Verify the repository name matches the configuration

### Styling Issues

1. Clear browser cache
2. Check that `custom.css` is properly referenced in `myst.yml`

## 📚 Resources

- [MyST Documentation](https://mystmd.org/)
- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [GitHub Pages Documentation](https://docs.github.com/en/pages)

## 👥 Contributors

- Joseph Loss
- Yuchen Duan
- Daniel Liberman

## 📄 License

This project is licensed under the MIT License.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

---

For questions or issues, please contact Joseph Loss at connect@josephjloss.com
