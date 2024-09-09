#!/usr/bin/env python3
# export_config.py

from config import config

def print_bash_exports(config):
    """Prints the configuration dictionary as bash export statements."""
    for key, value in config.items():
        # Convert Python lists to space-separated strings for bash arrays
        if isinstance(value, list):
            value = " ".join(map(str, value))
        print(f"export {key}='{value}'")

# Output the bash export statements
if __name__ == "__main__":
    print_bash_exports(config)

