#!/bin/bash

uv run myst clean --all

# Set the port for our local Jupyter process\np
port="9090"

# Define environment variables that will be used by MyST
# We'll use the values of these variables in our Jupyter server as well.
export JUPYTER_BASE_URL="http://localhost:${port}"
export JUPYTER_TOKEN="1234"

uv run jupyter-server --allow-root --ip 0.0.0.0 --port 9090 --IdentityProvider.token='1234' --ServerApp.allow_origin='*' &
echo $JUPYTER_BASE_URL

uv run myst build --execute --all --typst --md --html --site --pdf -d
#uv run myst start  --port 8265 --execute -d --server-port 9988
uv run myst start --execute -d


## jupytext backup
#uv run jupytext --from hyperparameter-optimization-whitepaper.ipynb --to hyperparameter-optimization-whitepaper.md  --execute
#uv run jupytext --from hyperparameter-optimization-whitepaper.ipynb --to md:myst,py:percent --sync --execute
#uv run jupytext --from hyperparameter-optimization-whitepaper.ipynb --to md:myst --sync --execute
#uv run jupytext --to md:myst --sync --execute hyperparameter-optimization-whitepaper.ipynb
#uv run jupytext --to md:myst --execute hyperparameter-optimization-whitepaper.ipynb
