"""Read the Hugo build used by the rendering checks."""

import os
from pathlib import Path

PUBLIC = Path(os.environ.get("GRASSMANNIAN_PUBLIC_DIR", Path(__file__).resolve().parents[1] / "public"))


def page(route):
    return (PUBLIC / route.strip("/") / "index.html").read_text()
