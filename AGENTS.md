# OpenCode Instructions

## Math and scientific notation

OpenCode uses KaTeX for mathematical rendering.

For all mathematical expressions:

- Inline math must use `$...$`.
- Display math must use `$$...$$`.
- Never use `\(...\)` or `\[...\]`.
- Prefer single-line display equations because multiline math rendering may be unreliable.
- Use standard KaTeX-compatible LaTeX syntax.
- Do not place Korean text inside math delimiters.
- Put explanatory prose outside the math environment.

Example:

The positional uncertainty is

$$\sigma_{\rm band}=\frac{{\rm FWHM}_{\rm band}}{2.355\,{\rm SNR}_{\rm min}}.$$

The matching criterion is

$$d<3\sqrt{\sigma_{\rm master}^2+\sigma_{\rm band}^2}.$$
