#!/usr/bin/env python3
"""Comment out entry(...) blocks in RMG training reactions files.

Usage:
  python comment_out_entries.py \
    --pattern "/home/moon/rmg/RMG-database/input/kinetics/families/Surface_*/training/reactions.py" \
    --backup

By default, modifies files in place. Use --dry-run to see planned changes.
"""

import argparse
import glob
import os


def iter_files(pattern):
    for path in glob.glob(pattern):
        if os.path.isfile(path):
            yield path


def comment_line(line):
    return "# " + line


def comment_entry_blocks(text):
    lines = text.splitlines(keepends=True)
    out = []
    in_block = False
    paren_balance = 0
    commented_blocks = 0

    for line in lines:
        if not in_block:
            if line.lstrip().startswith("entry(") or line.lstrip().startswith("entry ("):
                in_block = True
                paren_balance = line.count("(") - line.count(")")
                out.append(comment_line(line))
                if paren_balance <= 0:
                    in_block = False
                    commented_blocks += 1
                continue
            out.append(line)
            continue

        # inside block
        paren_balance += line.count("(") - line.count(")")
        out.append(comment_line(line))
        if paren_balance <= 0:
            in_block = False
            commented_blocks += 1

    return "".join(out), commented_blocks


def process_file(path, dry_run, backup):
    with open(path, "r", encoding="utf-8") as handle:
        original = handle.read()
    updated, blocks = comment_entry_blocks(original)
    changed = updated != original

    if changed and not dry_run:
        if backup:
            backup_path = path + ".bak"
            with open(backup_path, "w", encoding="utf-8") as handle:
                handle.write(original)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(updated)

    return blocks, changed


def main():
    parser = argparse.ArgumentParser(description="Comment out entry(...) blocks.")
    parser.add_argument(
        "--pattern",
        required=True,
        help="Glob pattern for reactions.py files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Do not modify files; only report changes.",
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Save a .bak copy before modifying.",
    )

    args = parser.parse_args()

    files = list(iter_files(args.pattern))
    if not files:
        print("No files matched pattern.")
        return 1

    total_blocks = 0
    changed_files = 0

    for path in files:
        blocks, changed = process_file(path, args.dry_run, args.backup)
        total_blocks += blocks
        changed_files += 1 if changed else 0
        status = "CHANGED" if changed else "UNCHANGED"
        print(f"{status}: {path} (commented {blocks} block(s))")

    print(f"Total files: {len(files)}; changed: {changed_files}; blocks commented: {total_blocks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
