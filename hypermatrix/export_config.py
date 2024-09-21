#!/usr/bin/env python3
# File: export_config.py

from config import config

def print_bash_exports(config):
    for key, value in config.items():
        if isinstance(value, list):
            value = " ".join(map(str, value))
        print(f"export {key}='{value}'")

if __name__ == "__main__":
    print_bash_exports(config)
