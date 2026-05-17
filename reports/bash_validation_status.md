# Bash Validation Status

- command_attempted: `bash --version`
- status: unavailable in this Windows environment
- observed_result: WSL/bash returned a Windows subsystem availability message rather than a usable shell.
- action: PowerShell quick-fixture validation was run and passed.
## Latest Check

- command: `bash --version`
- status: unavailable
- reason: WSL/bash is not configured in this Windows environment; PowerShell quick fixture validation passed instead.
