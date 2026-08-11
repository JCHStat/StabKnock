# Contributing

Contributions are welcome through GitHub issues and pull requests.

1. Do not commit restricted or participant-level data.
2. Keep scripts runnable from the repository root.
3. Use a fixed random seed for stochastic analyses.
4. Parse all R files before submitting: `Rscript -e "invisible(lapply(list.files('.', pattern='[.]R$', recursive=TRUE, full.names=TRUE), parse))"`.
5. Explain scientific or statistical changes in the pull request description.
