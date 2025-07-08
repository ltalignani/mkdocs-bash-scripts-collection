# cluster-usage

## mkdocs-material Initial Installation

1. Open a terminal in the folder that you want to create the project in
2. Type which python or which python3 to check Python is installed
3. Setup a Python virtual environment by typing `python -m venv venv`
4. Type `source venv/bin/activate` to activate the virtual environment
5. Install mkdocs material - `pip install mkdocs-material`
6. Open Visual studio code in this folder with `code .`
7. Open a terminal within Visual Studio code
8. Activate the virtual environment on the new terminal with `source venv/bin/activate`
9. Create a new site `mkdocs new .`
10. Add the following basic mkdocs.yml configuration:
```yaml
site_name: My MkDocs Material Documentation
site_url: https://sitename.example
theme:
  name: material
```

11. Type `mkdocs serve` to launch the site
12. Check the site at http://localhost:8000