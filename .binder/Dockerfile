FROM jupyter/base-notebook:python-3.11

ENV ZSH="$HOME/.oh-my-zsh" \
		PYTHONUNBUFFERED=1 \
		PYTHONDONTWRITEBYTECODE=1 \
 		PIP_NO_CACHE_DIR=OFF \
    PIP_DEFAULT_TIMEOUT=100 \
		PATH="$PATH:$HOME/bin/" \
		UV_VENV=.venv \
		UV_USE_PYPACKAGES=1 \
    MSYS_NO_PATHCONV=1 

#export ZSH_AUTOSUGGEST_STRATEGY="history"
#fpath+=${ZSH_CUSTOM:-${ZSH:-~/.oh-my-zsh}/custom}/plugins/zsh-completions/src
#source $ZSH/oh-my-zsh.sh

# ENV NVM_DIR="$HOME/.nvm"
# [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"
# [ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"

# alias for env
#. "$HOME/.local/bin/env"

# Switch to root for system installations
USER root

# Install system dependencies from apt.txt and additional required packages
RUN apt-get update && apt-get install -y \
    ack \
    alien \
    apt-utils \
    bat \
    build-essential \
    cargo \
    cmake \
    curl \
    fzf \
    fonts-liberation \
    gcc \
    graphviz \
    imagemagick \
    inkscape \
    jq \
    latexmk \
    libatk-bridge2.0-0 \
    libgtk-3-0 \
    libssl-dev \
    libx11-xcb1 \
    libxcomposite1 \
    libxdamage1 \
    libxrandr2 \
    nano \
    nodejs \
    npm \
    pandoc \
    pkg-config \
    python3-dev \
    ripgrep \
    rpm \
    texinfo \
    texlive-fonts-recommended \
    texlive-plain-generic \
    texlive-xetex \
    tldr \
    webp \
    wkhtmltopdf \
    xclip \
    xdg-utils \
    zsh \
    && rm -rf /var/lib/apt/lists/*

# Install Rust using rustup (latest stable version)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && \
    echo 'source "$HOME/.cargo/env"' >> ~/.bashrc && \
    export PATH="$HOME/.cargo/bin:$PATH"
ENV PATH="/root/.cargo/bin:${PATH}"

# Copy all Rust binaries to /usr/local/bin to make available for all users
RUN . "$HOME/.cargo/env" && \
    cp -a $HOME/.cargo/bin/* /usr/local/bin/ && \
    chmod 755 /usr/local/bin/cargo

# Install Oh My Zsh, zsh plugins, and set up zsh as default shell for jovyan user
RUN apt-get update && apt-get install -y --no-install-recommends \
    fzf \
    bat \
    jq \
    tldr \
    nano \
    && ln -s /usr/bin/batcat /usr/local/bin/bat \
    && rm -rf /var/lib/apt/lists/*

# Install Oh My Zsh for jovyan
USER jovyan
ENV ZSH="$HOME/.oh-my-zsh"
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)" "" --unattended \
    && git clone https://github.com/zsh-users/zsh-autosuggestions ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-autosuggestions \
    && git clone https://github.com/zsh-users/zsh-completions ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-completions \
    && mkdir -p ~/bin ~/.local/bin

# Set up zsh as default shell for jovyan user
USER root
RUN chsh -s /bin/zsh jovyan

# Install Typst - note that it installs to /home/jovyan/.typst
RUN curl -fsSL https://typst.community/typst-install/install.sh | sh && \
    # Since it's installed to jovyan's home directory, make it accessible to all
    mkdir -p /usr/local/typst/bin && \
    cp /home/jovyan/.typst/bin/typst /usr/local/typst/bin/ && \
    chmod -R 755 /usr/local/typst && \
    # Add to system PATH
    echo 'export PATH="/usr/local/typst/bin:$PATH"' >> /etc/profile && \
    echo 'export PATH="/usr/local/typst/bin:$PATH"' >> /etc/bash.bashrc

# Install cargo packages system-wide and ensure they're in PATH
RUN . "$HOME/.cargo/env" && \
    cargo install --locked typst-cli gitui eza \
    && cp ~/.cargo/bin/* /usr/local/bin/

# Corrected command to install uv and uvx
# 1. Download and run the uv installation script. It installs to ~/.local/bin by default.
# 2. Move uv from ~/.local/bin/uv to /usr/local/bin/uv.
# 3. Move uvx from ~/.local/bin/uvx to /usr/local/bin/uvx.
# 4. Set execute permissions for both.
RUN curl -fsSL https://astral.sh/uv/install.sh | sh && \
    mv ~/.local/bin/uv /usr/local/bin/uv && \
    mv ~/.local/bin/uvx /usr/local/bin/uvx && \
    chmod 755 /usr/local/bin/uv && \
    chmod 755 /usr/local/bin/uvx

# Set environment variables for Puppeteer
ENV PUPPETEER_EXECUTABLE_PATH="/usr/bin/chromium"
ENV PATH="/usr/local/typst/bin:${PATH}"

# Copy the repository content to the home directory
# This ensures proper ownership of files
COPY --chown=${NB_UID}:${NB_GID} . ${HOME}

# Switch to the repository directory
WORKDIR ${HOME}

# Install npm packages in the repository directory
# This ensures build-components.mjs can find them
RUN npm install --force

# Switch back to jovyan user for remaining operations
USER ${NB_UID}

# Install Python packages using uv
WORKDIR ${HOME}
RUN /usr/local/bin/uv pip install --system --no-cache-dir -r .binder/requirements.txt

# Install NVM for Node.js version management
RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.7/install.sh | bash

#todo: Create env file with any additional environment variables
# RUN mkdir -p ${HOME}/.local/bin && \
#     echo '#!/bin/bash' > ${HOME}/.local/bin/env && \
#     echo 'echo "Environment loaded"' >> ${HOME}/.local/bin/env && \
#     chmod +x ${HOME}/.local/bin/env

# Configure jovyan's zshrc with the exact requested settings
#COPY .zshrc ${HOME}/.zshrc

# Ensure tools are available in all shells by adding them to global profile
# RUN echo 'export PATH="/usr/local/bin:/usr/local/typst/bin:$HOME/.local/bin:$HOME/.cargo/bin:$HOME:$PATH"' >> /etc/profile && \
#     echo 'export PATH="/usr/local/bin:/usr/local/typst/bin:$HOME/.local/bin:$HOME/.cargo/bin:$HOME:$PATH"' >> /etc/bash.bashrc && \
#     echo 'source "$HOME/.cargo/env"' >> /etc/profile && \
#     echo 'source "$HOME/.cargo/env"' >> /etc/bash.bashrc

# Build components using Node.js
RUN node build-components.mjs || echo "build-components.mjs failed"

# Build JupyterLab using uv
RUN /usr/local/bin/uv run --no-cache jupyter lab build --dev-build=False || echo "jupyter lab build failed"

# Set the working directory to home
WORKDIR ${HOME}

# Ensure all tools are in path when container starts
#RUN echo 'source /etc/profile' >> /etc/zsh/zshenv

# Command to run JupyterLab using uv
CMD ["sh", "-c", "/usr/local/bin/uv run jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root --NotebookApp.token='' --NotebookApp.password=''"]
