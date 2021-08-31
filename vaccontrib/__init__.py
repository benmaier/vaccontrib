# -*- coding: utf-8 -*-
"""
Initializes this package with metadata.
"""

from .metadata import (
        __version__,
        __author__,
        __copyright__,
        __credits__,
        __license__,
        __maintainer__,
        __email__,
        __status__,
    )

from .main import (
        get_next_generation_matrix_from_matrices,
        get_contribution_matrix,
        get_reduced_contribution_matrix,
        get_reduced_vaccinated_susceptile_contribution_matrix,
        get_reduced_population_contribution_matrix,
    )
