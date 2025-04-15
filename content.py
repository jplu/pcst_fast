from pathlib import Path

sources_paths = [p for p in Path("src").rglob("**/*") if p.is_file() and p.as_posix().endswith(".cc")]
includes_paths = [p for p in Path("include").rglob("**/*") if p.is_file() and p.as_posix().endswith(".h")]
bindings_paths = [p for p in Path("bindings").rglob("**/*") if p.is_file() and p.as_posix().endswith(".cc")]
final_content = ""

for p in sources_paths:
    final_content += f"`{p}`:\n"
    final_content += f"```\n{p.read_text()}\n```\n"

for p in includes_paths:
    final_content += f"`{p}`:\n"
    final_content += f"```\n{p.read_text()}\n```\n"

for p in sources_paths:
    final_content += f"`{p}`:\n"
    final_content += f"```\n{p.read_text()}\n```\n"

final_content += f"`setup.py`:\n```\n{Path('setup.py').read_text()}\n```\n"
final_content += f"`Makefile`:\n```\n{Path('Makefile').read_text()}\n```\n"
final_content += f"`pyproject.toml`:\n```\n{Path('pyproject.toml').read_text()}\n```\n"

Path("content.txt").write_text(final_content)