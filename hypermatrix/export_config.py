#!/usr/bin/env python3
# export_config.py

from config import config
import logging, logging.handlers # add robust logging to all files



# import logging, logging.handlers
# cfg.add_namespace(logging)
# print cfg.logging.loggers['input.xls'].handlers[0][0]
# print cfg.logging.handlers.console[1].stream
# print cfg['logging']['handlers']['console'][1]['stream']
# print cfg.logging.handlers.email[1]['from']
# x = cfg.logging.handlers.email[1].to
# print x
# for a in x:
#     print a
# print x[0:2]
# print cfg.logging.handlers.file[1].filename



print("source: export_config.py")
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
