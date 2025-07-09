# Options Pricing in Lévy Models - MyST Site

This repository contains a MyST (Markedly Structured Text) multi-page documentation site for the Options Pricing in Lévy Models research project, properly configured with the book theme.

## 📁 Project Structure

```
option-pricing-myst-site/
├── .github/
│   └── workflows/
│       └── deploy-myst.yml          # GitHub Actions workflow for deployment
├── _static/
│   └── custom.css                   # Custom styling for book theme
├── assets/                          # Directory for images
├── documentation/                   # PDF documents and slides
│   ├── misc/                        # Slide images (Slide1.PNG - Slide7.PNG)
│   ├── Project Research Proposal.pdf
│   ├── Option Pricing in Levy Models - Feng et al - Academic Paper.pdf
│   └── R-Finance Presentation Slides.pdf
├── exports/                         # Export configurations
├── .gitignore                       # Git ignore file
├── DEPLOYMENT_CHECKLIST.md          # Step-by-step deployment guide
├── README.md                        # This file
├── myst.yml                         # MyST configuration (book theme)
├── index.md                         # Page 1: Research Proposal
├── academic-paper.md                # Page 2: Academic Paper
└── presentation-slides.md           # Page 3: R/Finance Presentation
```

## 🚀 Deployment Instructions

### Prerequisites

1. **GitHub Repository**: Ensure your code is pushed to GitHub
2. **GitHub Pages**: Enable GitHub Pages in your repository settings
3. **PDFs and Images**: Place your files in the correct directories

### Step 1: Enable GitHub Pages

1. Go to your repository on GitHub
2. Navigate to **Settings** → **Pages**
3. Under **Source**, select **GitHub Actions**

### Step 2: Add Required Files

Place these files in the `documentation/` directory:
- `Project Research Proposal.pdf`
- `Option Pricing in Levy Models - Feng et al - Academic Paper.pdf`
- `R-Finance Presentation Slides.pdf`

Place slide images in `documentation/misc/`:
- `Slide1.PNG` through `Slide7.PNG`

### Step 3: Push to GitHub

```bash
# Initialize git (if not already done)
git init

# Add all files
git add .

# Commit
git commit -m "Initial MyST book theme site for Options Pricing research"

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
2. Update `myst.yml` if adding new pages to the table of contents
3. Push changes to trigger automatic deployment

### MyST Syntax Features

The site uses MyST markdown with book theme features:

- **Admonitions**: 
  ```markdown
  :::{admonition} Note
  :class: note
  
  Your content here
  :::
  ```

- **Full-width content**:
  ```markdown
  :::{.full-width}
  Your wide content here
  :::
  ```

- **Cards**:
  ```markdown
  :::{card}
  **Title**
  ^^^
  Content
  :::
  ```

- **Button links**:
  ```markdown
  ```{button-link} path/to/file.pdf
  :color: primary
  :align: center
  📄 Download PDF
  ```

### Adding Images

1. Place images in the appropriate directory
2. Reference in markdown using MyST figure directive:
   ```markdown
   ```{figure} path/to/image.png
   :width: 100%
   :align: center
   :alt: Image description
   
   Caption for the image
   ```

## 🎨 Book Theme Features

The site uses MyST's book theme with:

- **Clean typography**: Optimized for readability
- **Responsive layout**: Works on all devices
- **Navigation sidebar**: Easy access to all pages
- **Professional styling**: Academic paper appearance
- **Code highlighting**: R code examples with syntax highlighting
- **Mathematical equations**: Full LaTeX support

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

- **Book theme layout**: Professional academic appearance
- **Multi-page structure**: Three interconnected pages with navigation
- **PDF integration**: Direct links to research documents
- **Slide presentation**: All 7 slides embedded in the presentation page
- **Responsive design**: Optimized for all devices
- **Code highlighting**: R code examples with line numbers
- **Mathematical equations**: LaTeX support for formulas
- **Interactive elements**: Cards, admonitions, and button links

## 🐛 Troubleshooting

### Build Failures

1. Check the Actions tab for error messages
2. Ensure all markdown files are valid MyST syntax
3. Verify PDF and image files exist in the correct locations

### Missing Images

1. Ensure slide images are named exactly: `Slide1.PNG` through `Slide7.PNG`
2. Place them in `documentation/misc/` directory
3. Check case sensitivity (PNG not png)

### Styling Issues

1. Clear browser cache
2. Verify `_static/custom.css` exists
3. Check that book theme is properly configured in `myst.yml`

## 📚 Resources

- [MyST Documentation](https://mystmd.org/)
- [MyST Book Theme](https://mystmd.org/guide/themes)
- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [GitHub Pages Documentation](https://docs.github.com/en/pages)

## 👥 Contributors

- Joseph Loss - Lead Researcher
- Yuchen Duan - Co-Researcher
- Daniel Liberman - Co-Researcher

Illinois Institute of Technology

## 📄 License

This project is licensed under the MIT License.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

---

For questions or issues, please contact Joseph Loss at connect@josephjloss.com
