# Optional Method Installation Status

This file records unavailable optional methods and does not claim successful execution.

- scvi-tools: not installed. Suggested command: `pip install scvi-tools`.
- TIGON: not installed. Suggested command: `pip install tigon`.
- TrajectoryNet: not installed. Suggested command: `pip install git+https://github.com/KrishnaswamyLab/TrajectoryNet.git`.
- MIOFlow: not installed. Suggested command: `pip install git+https://github.com/KrishnaswamyLab/MIOFlow.git`.
- GEARS: not installed. Suggested command: `pip install gears`.
- scGPT: not installed. Suggested command: `pip install scgpt`.
- scFoundation: not installed. Suggested command: `install from the maintained scFoundation repository/release`.

Dry-run attempts recorded during this build:

1. Sandboxed `python -m pip install --dry-run scvi-tools tigon gears mioflow` failed with Windows socket permission/network sandbox errors.
2. Escalated non-mutating dry-run reached PyPI, downloaded scvi-tools metadata, then failed because `tigon` is not available on PyPI as `tigon`.

