"""Package the five storyboard figures, evidence and regeneration source."""
from pathlib import Path
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/eight_module_storyboards"
OUTPUT = ROOT / "output/figures/eight_module_storyboards.zip"
manifest = json.loads((RESULTS / "manifest.json").read_text())
entries = {}

for path in sorted(RESULTS.rglob("*")):
    if path.is_file():
        entries[str(path.relative_to(RESULTS))] = path.read_bytes()

captions = "\n\n".join(
    f"FIGURE {i}: {study['title']}\n\n{study['caption']}"
    for i, study in enumerate(manifest["studies"], 1)
)
entries["captions.txt"] = (captions + "\n").encode()
entries["speaker_script.txt"] = (
    ROOT / "output/presentations/electroexocytosis_eight_module_study_scripts.txt"
).read_bytes()

source_paths = [ROOT / "pyproject.toml", ROOT / "README.md"]
source_paths += [
    path for path in (ROOT / "electro_exocytosis").rglob("*")
    if path.is_file()
    and path.suffix in {".py", ".yaml", ".yml", ".csv"}
    and "assets" not in path.relative_to(ROOT).parts
]
source_paths += [
    ROOT / "examples/__init__.py",
    ROOT / "examples/storyboard_eight_modules.py",
    ROOT / "examples/storyboard_multilayer_overview.py",
    ROOT / "tests/test_eight_module_storyboards.py",
]
for path in source_paths:
    if path.exists():
        entries["code/" + str(path.relative_to(ROOT))] = path.read_bytes()

checksums = "".join(
    f"{hashlib.sha256(data).hexdigest()}  {name}\n"
    for name, data in sorted(entries.items())
)
entries["SHA256SUMS.txt"] = checksums.encode()
with zipfile.ZipFile(OUTPUT, "w", zipfile.ZIP_DEFLATED, compresslevel=7) as archive:
    for name, data in sorted(entries.items()):
        archive.writestr("eight_module_storyboards/" + name, data)
with zipfile.ZipFile(OUTPUT) as archive:
    assert archive.testzip() is None
    for name, data in entries.items():
        assert archive.read("eight_module_storyboards/" + name) == data
print(f"Packaged {len(entries)} files in {OUTPUT} ({OUTPUT.stat().st_size:,} bytes)")
