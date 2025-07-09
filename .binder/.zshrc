ZSH_THEME="robbyrussell"
plugins=(
  git
  zsh-autosuggestions
#  F-Sy-H
  tldr
)

#export ZSH_AUTOSUGGEST_STRATEGY="history"
fpath+=${ZSH_CUSTOM:-${ZSH:-~/.oh-my-zsh}/custom}/plugins/zsh-completions/src
source $ZSH/oh-my-zsh.sh

export PATH="$PATH:$HOME/bin/"
export UV_VENV=.venv
export UV_USE_PYPACKAGES=1
export MSYS_NO_PATHCONV=1
export PUPPETEER_EXECUTABLE_PATH="/usr/bin/chromium"

# alias for env
. "$HOME/.local/bin/env"
alias j!=jbang
alias ls="ls --group-directories-first -a -h -p --format=vertical -s --color=always"
alias bat="batcat"
alias fcat="fzf --preview \"batcat --style=numbers --color=always --line-range=:500 {}\""
alias zshconfig="nano ~/.zshrc"
alias cht="cht.sh"

export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"
[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"

export PATH="$HOME/.jbang/bin:$HOME/.jbang/currentjdk/bin:/usr/local/bin:$PATH"
export JAVA_HOME=$HOME/.jbang/currentjdk
