#!/usr/bin/env python3
"""Copy a powheg.input, replacing the VALUE of named keys.

    set_keys.py IN OUT key=value [key=value ...]

Only the value token after the key at the start of a line is replaced; the
trailing "! comment" is kept.  A key that is not present in the card is an
ERROR, not an append: POWHEG reads an absent key as its own default without a
word (CONVENTIONS.md rule 2), so a typo such as `ih2=` against a card spelling it
`ih2 ` must stop the driver rather than silently run the proton.

IN and OUT may be the same file.  OUT is written as a real file, never
through a symlink (see powheg_install_card in powheg/powheg_variants.sh).
"""
import os
import re
import sys


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    text = open(src).read()
    for kv in sys.argv[3:]:
        key, _, val = kv.partition("=")
        if not key or not val:
            raise SystemExit(f"set_keys.py: bad assignment {kv!r}")
        pat = re.compile(r"^(" + re.escape(key) + r")([ \t]+)(\S+)", re.M)
        text, n = pat.subn(lambda m: m.group(1) + m.group(2) + val, text)
        if n != 1:
            raise SystemExit(f"set_keys.py: key {key!r} found {n} times in "
                             f"{src} (want exactly 1)")
    if os.path.islink(dst):
        os.remove(dst)
    with open(dst, "w") as f:
        f.write(text)


if __name__ == "__main__":
    main()
