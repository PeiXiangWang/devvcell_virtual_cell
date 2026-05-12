from __future__ import annotations

import torch


def barycentric_loss(pred_delta: torch.Tensor, target_delta: torch.Tensor) -> torch.Tensor:
    return torch.mean((pred_delta - target_delta) ** 2)


def smoothness_loss(pred_delta: torch.Tensor) -> torch.Tensor:
    return 1e-4 * torch.mean(pred_delta**2)

