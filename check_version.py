#!/usr/bin/env python

import re
import requests


pypi_version = requests.get("https://test.pypi.org/pypi/atlas4py/json").json()["info"]["version"]

with open("src/atlas4py/_version.py", "r") as file:
    match = re.search("""__version__[ ]*=[ ]*[\'\"]([^ \'\"]+)""", file.read().replace("\n", " "))
    local_version = match[1].strip()


num_re = re.compile("[a-zA-Z]*([0-9]+)")
local_parts = [int(num_re.match(p)[1]) for p in local_version.split(".")]
pypi_parts = [int(num_re.match(p)[1]) for p in pypi_version.split(".")]

if local_parts == pypi_parts or not all(l >= p for l, p in zip(local_parts, pypi_parts)):
    raise ValueError(
        f"Local version '{local_version}'' is not newer than the currently published version '{pypi_version}'"
    )

print("Ok")
