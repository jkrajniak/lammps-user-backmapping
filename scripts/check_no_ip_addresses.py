#!/usr/bin/env python3
"""Pre-commit guard: reject literal IPv4 addresses in tracked files.

Infrastructure IPs (VM hosts, etc.) don't belong in source control -- they
leak internal topology and go stale the moment a VM is recreated. Point at
a placeholder (`<vm_host>`) or an env var / config file instead.
"""

from __future__ import annotations

import re
import sys

IPV4_RE = re.compile(r"\b(?:\d{1,3}\.){3}\d{1,3}\b")

# Not "an IP address of a VM" in the sense this hook cares about.
ALLOWLIST = {
    "0.0.0.0",
    "127.0.0.1",
    "255.255.255.255",
    "8.8.8.8",  # common DNS example
}


def octets_valid(ip: str) -> bool:
    return all(0 <= int(part) <= 255 for part in ip.split("."))


def find_ips(text: str) -> list[str]:
    hits = []
    for match in IPV4_RE.finditer(text):
        ip = match.group(0)
        if ip in ALLOWLIST or not octets_valid(ip):
            continue
        hits.append(ip)
    return hits


def main(argv: list[str]) -> int:
    failed = False
    for path in argv:
        try:
            with open(path, encoding="utf-8", errors="ignore") as f:
                text = f.read()
        except OSError:
            continue
        ips = find_ips(text)
        if ips:
            failed = True
            unique = sorted(set(ips))
            print(f"{path}: literal IP address(es) found: {', '.join(unique)}")
    if failed:
        print(
            "\nDon't commit literal IP addresses (VM hosts, etc.) -- use a "
            "placeholder like <vm_host> or read it from config/env instead."
        )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
