from pathlib import Path
import sys
import re
import sysconfig


TARGET_FILES = [
    Path("/opt/venv/lib/python3.11/site-packages/scanet/preprocessing.py"),
    Path("/opt/venv/lib/python3.11/site-packages/scanet/coexpression.py"),
    Path("/opt/venv/lib/python3.11/site-packages/scanet/plottingNetworks.py")
]

PATTERNS_TO_REMOVE = [
    "from turtle import up\n",
    "from turtle import up\r\n",
    "import turtle\n",
    "import turtle\r\n",
    "from turtle import title\n",
    "from turtle import title\r\n",
    "from turtle import shape\n",
    "from turtle import shape\r\n",
]


def patch_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"[setup_docker] File not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    for pattern in PATTERNS_TO_REMOVE:
        text = text.replace(pattern, "")

    if text != original:
        path.write_text(text, encoding="utf-8")
        print(f"[setup_docker] Patched {path}: removed turtle import(s)")
    else:
        print(f"[setup_docker] No turtle import found in {path}")


def main() -> None:
    for path in TARGET_FILES:
        patch_file(path)
    # PySCENIC 0.11 uses NumPy aliases removed in NumPy 1.24.
    # Replace only those exact deprecated spellings in the installed package.
    package = Path(sysconfig.get_paths()["purelib"]) / "pyscenic"
    for path in package.rglob("*.py"):
        text = path.read_text()
        updated = re.sub(r"\bnp\.(object|float|int|bool|str|complex)\b", r"\1", text)
        if updated != text:
            path.write_text(updated)


if __name__ == "__main__":
    main()
