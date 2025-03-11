# How to write this paper

## General process

1. Edit the file `main.md` (You can do this in the web editor!)
2. Commit your changes to the main branch (You can do this in the web editor!)

That's it!

GitHub actions will take care of updating the `.tex`, `.docx` and `.pdf` files.

You can also open a pull request, for example by selecting "Create a new branch for this commit and start a pull request". This provides a nice interface for discussing changes to the manuscript with co-authors.

## How to preview changes locally

To preview your edits locally, e.g. in order to check if your table or equation renders correctly, you will need to install [pandoc](https://pandoc.org/MANUAL.html) and a latex engine (e.g. [tectonic](https://tectonic-typesetting.github.io/en-US/)).

Then go into this folder (i.e. `cmfapoc/paper`) and run a command like this:

```shell
pandoc --standalone -o main.pdf --template=arxiv-template.tex --pdf-engine=tectonic --natbib main.md
```
