#!/bin/bash
find /home/agent -name "geneformer" -type d 2>/dev/null
find /usr -path "*/geneformer/__init__.py" 2>/dev/null  
find /home/agent -path "*/geneformer/__init__.py" 2>/dev/null
ls /home/agent/.local/lib/python*/site-packages/geneformer/ 2>/dev/null
python3.12 -c "import geneformer; print(geneformer.__file__)" 2>/dev/null
