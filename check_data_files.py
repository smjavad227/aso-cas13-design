import os
from pathlib import Path

data_dir = Path("outputs/data")
if not data_dir.exists():
    print("Directory outputs/data does not exist")
else:
    for f in data_dir.iterdir():
        if f.is_file():
            size = f.stat().st_size
            with f.open() as fh:
                try:
                    first = next(fh).strip()
                except StopIteration:
                    first = "<EMPTY>"
            print(f"{f.name}\t{size} bytes\tFirst line: {first}")
