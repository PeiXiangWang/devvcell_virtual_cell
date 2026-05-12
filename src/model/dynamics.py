from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


class PositiveScalar(torch.nn.Module):
    def __init__(self, init: float):
        super().__init__()
        self.raw = torch.nn.Parameter(torch.tensor(float(init)).log().expm1())

    def forward(self) -> torch.Tensor:
        return F.softplus(self.raw)


class SmallMLP(torch.nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, output_dim: int):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(input_dim, hidden_dim),
            torch.nn.SiLU(),
            torch.nn.Linear(hidden_dim, hidden_dim),
            torch.nn.SiLU(),
            torch.nn.Linear(hidden_dim, output_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


@dataclass(frozen=True)
class DynamicsFlags:
    use_teacher: bool = False
    use_swarm: bool = False
    use_birth_death: bool = False
    use_diffusion: bool = False
    use_cci: bool = False
    use_memory: bool = False


class SwarmLineageDynamics(torch.nn.Module):
    """Trainable finite-agent dynamics used by SwarmLineage-OT variants."""

    def __init__(self, feature_dim: int, latent_dim: int, hidden_dim: int = 96):
        super().__init__()
        self.latent_dim = int(latent_dim)
        self.v_intrinsic_net = SmallMLP(feature_dim, hidden_dim, latent_dim)
        self.v_teacher_net = SmallMLP(feature_dim, hidden_dim, latent_dim)
        self.birth_net = SmallMLP(feature_dim, hidden_dim // 2, 1)
        self.death_net = SmallMLP(feature_dim, hidden_dim // 2, 1)
        self.sigma_net = SmallMLP(feature_dim, hidden_dim // 2, 1)
        self.cci_gate = SmallMLP(feature_dim, hidden_dim // 2, 1)
        self.alignment_weight = PositiveScalar(0.08)
        self.cohesion_weight = PositiveScalar(0.05)
        self.separation_weight = PositiveScalar(0.03)
        self.memory_weight = PositiveScalar(0.03)
        self.max_step = 3.0

    def vector_field(
        self,
        features: torch.Tensor,
        swarm_delta: torch.Tensor | None = None,
        cci_delta: torch.Tensor | None = None,
        memory_delta: torch.Tensor | None = None,
        flags: DynamicsFlags = DynamicsFlags(),
    ) -> torch.Tensor:
        v = self.v_intrinsic_net(features)
        if flags.use_teacher:
            v = v + self.v_teacher_net(features)
        if flags.use_swarm and swarm_delta is not None:
            v = v + (self.alignment_weight() + self.cohesion_weight() + self.separation_weight()) * swarm_delta
        if flags.use_cci and cci_delta is not None:
            gate = torch.sigmoid(self.cci_gate(features))
            v = v + gate * cci_delta
        if flags.use_memory and memory_delta is not None:
            v = v + self.memory_weight() * memory_delta
        return torch.clamp(v, min=-self.max_step, max=self.max_step)

    def birth_hazard(self, features: torch.Tensor, flags: DynamicsFlags) -> torch.Tensor:
        if not flags.use_birth_death:
            return torch.zeros(features.shape[0], device=features.device)
        return F.softplus(self.birth_net(features)).squeeze(-1)

    def death_hazard(self, features: torch.Tensor, flags: DynamicsFlags) -> torch.Tensor:
        if not flags.use_birth_death:
            return torch.zeros(features.shape[0], device=features.device)
        return F.softplus(self.death_net(features)).squeeze(-1)

    def sigma(self, features: torch.Tensor, flags: DynamicsFlags) -> torch.Tensor:
        if not flags.use_diffusion:
            return torch.full((features.shape[0],), 0.015, device=features.device)
        return 0.005 + 0.18 * torch.sigmoid(self.sigma_net(features)).squeeze(-1)

