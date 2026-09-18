from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.inset import InsetIndicator
from matplotlib.text import Text
from PIL import Image
from pypdf import PdfReader
import pytest


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "scripts" / "generate_jcp_figures_v34.py"


def _load_generator():
    spec = importlib.util.spec_from_file_location("generate_jcp_figures_v34", GENERATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _capture_generated_figure(generator, monkeypatch, tmp_path: Path, function_name: str):
    """Run the real figure function while intercepting only filesystem output."""
    captured = {}

    def capture_figure(fig, stem: str) -> None:
        captured[stem] = fig

    monkeypatch.setattr(generator, "save", capture_figure)
    monkeypatch.setattr(generator, "DATA", tmp_path)
    getattr(generator, function_name)()
    assert len(captured) == 1
    return next(iter(captured.values()))


def _visible_text_sizes(fig) -> list[float]:
    fig.canvas.draw()
    return [
        float(text.get_fontsize())
        for text in fig.findobj(match=Text)
        if text.get_visible() and text.get_text().strip()
    ]


def test_emitted_file_metadata_is_publication_facing(
    monkeypatch, tmp_path: Path
) -> None:
    """Upload-facing figure files must not expose stale baseline branding."""
    generator = _load_generator()
    monkeypatch.setattr(generator, "FIG", tmp_path)

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    generator.save(fig, "metadata_probe")

    pdf_metadata = PdfReader(tmp_path / "metadata_probe.pdf").metadata
    assert pdf_metadata is not None
    assert pdf_metadata["/Title"] == (
        "Complete-positivity topology of Runge-Kutta maps: publication figures"
    )
    assert pdf_metadata["/Creator"] == (
        "JCP deterministic Matplotlib figure generator"
    )
    assert pdf_metadata["/CreationDate"] == "D:20260820000000Z"
    assert pdf_metadata["/ModDate"] == "D:20260820000000Z"

    with Image.open(tmp_path / "metadata_probe.png") as image:
        assert image.info["Software"] == (
            "JCP deterministic Matplotlib figure generator"
        )
        assert image.info["Creation Time"] == "2026-08-20T00:00:00Z"


def test_buffered_component_legend_does_not_cover_any_component_bar(
    monkeypatch, tmp_path: Path
) -> None:
    """The panel-a legend must remain readable without hiding plotted intervals."""
    generator = _load_generator()
    captured = {}

    def capture_figure(fig, stem: str) -> None:
        captured[stem] = fig

    monkeypatch.setattr(generator, "save", capture_figure)
    monkeypatch.setattr(generator, "DATA", tmp_path)
    generator.fig2_robustness_noncommuting()

    fig = captured["Fig2_robustness"]
    try:
        fig.canvas.draw()
        ax = fig.axes[0]
        renderer = fig.canvas.get_renderer()
        legend = ax.get_legend()
        assert legend is not None
        legend_box = legend.get_window_extent(renderer)

        covered = []
        for component_bar in ax.collections:
            data_box = component_bar.get_datalim(ax.transData)
            display_box = ax.transData.transform_bbox(data_box).padded(5)
            if legend_box.overlaps(display_box):
                covered.append(data_box.bounds)

        assert not covered, f"legend overlaps component bars: {covered}"
    finally:
        plt.close(fig)


def test_main_figures_are_authored_at_journal_width_with_readable_type(
    monkeypatch, tmp_path: Path
) -> None:
    """Figures 1--5 must not rely on severe manuscript down-scaling for legibility."""
    generator = _load_generator()
    functions = (
        "fig1_topology_conditioning",
        "fig2_robustness_noncommuting",
        "fig3_candidate_maps",
        "fig4_certified_guard_conditioning",
        "fig5_tolerance_sweep",
    )

    for function_name in functions:
        fig = _capture_generated_figure(
            generator, monkeypatch, tmp_path, function_name
        )
        try:
            width_inches, _ = fig.get_size_inches()
            assert width_inches <= 7.6, (
                f"{function_name} is {width_inches:.2f} in wide; "
                "finished-size text would be reduced below journal guidance"
            )
            sizes = _visible_text_sizes(fig)
            assert sizes and min(sizes) >= 7.5, (
                f"{function_name} contains visible text at {min(sizes):.2f} pt"
            )
        finally:
            plt.close(fig)


def test_figure5_equal_work_pair_is_explained_and_series_are_shape_coded(
    monkeypatch, tmp_path: Path
) -> None:
    """The equal-work RK4 pair must be explicit and all benchmark series distinguishable without color."""
    generator = _load_generator()
    fig = _capture_generated_figure(
        generator, monkeypatch, tmp_path, "fig5_tolerance_sweep"
    )
    try:
        accuracy_ax, action_ax = fig.axes[:2]
        labels = {
            "error-only RK4",
            "CP guard + rotating fallback",
            "rotating-frame RK4",
        }
        series = [line for line in accuracy_ax.lines if line.get_label() in labels]
        assert {line.get_label() for line in series} == labels
        assert len({line.get_marker() for line in series}) == 3
        assert len({line.get_linestyle() for line in series}) == 3

        explanatory_text = " ".join(
            text.get_text() for text in accuracy_ax.texts
        ).replace(" ", "")
        assert "10^{-2}" in explanatory_text
        assert "3\\times10^{-3}" in explanatory_text

        hatches = {
            patch.get_hatch()
            for container in action_ax.containers
            for patch in container.patches
        }
        assert len(hatches) == 3
        assert all(hatches)
    finally:
        plt.close(fig)


def test_figure5_roundoff_reference_does_not_compress_accuracy_curves(
    monkeypatch, tmp_path: Path
) -> None:
    """A roundoff-floor comparator must not compress or clutter the log panel."""
    generator = _load_generator()
    fig = _capture_generated_figure(
        generator, monkeypatch, tmp_path, "fig5_tolerance_sweep"
    )
    try:
        accuracy_ax = fig.axes[0]
        legend_labels = accuracy_ax.get_legend_handles_labels()[1]
        assert "Strang (commuting reference)" not in legend_labels

        panel_text = " ".join(text.get_text() for text in accuracy_ax.texts).lower()
        assert "strang" not in panel_text

        lower, upper = accuracy_ax.get_ylim()
        assert lower >= 1e-7
        assert upper / lower <= 1e5
    finally:
        plt.close(fig)


def test_figure5_omits_all_zero_exact_fallback_stack(
    monkeypatch, tmp_path: Path
) -> None:
    """A candidate type with zero accepts in every scenario must not occupy the legend."""
    generator = _load_generator()
    fig = _capture_generated_figure(
        generator, monkeypatch, tmp_path, "fig5_tolerance_sweep"
    )
    try:
        action_ax = fig.axes[1]
        legend_labels = action_ax.get_legend_handles_labels()[1]
        assert legend_labels == ["direct", "projected", "rotating fallback"]
        assert len(action_ax.containers) == 3

        panel_text = " ".join(text.get_text() for text in action_ax.texts).lower()
        assert "exact fallback" not in panel_text
    finally:
        plt.close(fig)


def test_figure5_annotations_clear_legend_curves_and_bars(
    monkeypatch, tmp_path: Path
) -> None:
    """Rendered annotations must not collide with the legend, blue curve, or bars."""
    generator = _load_generator()
    fig = _capture_generated_figure(
        generator, monkeypatch, tmp_path, "fig5_tolerance_sweep"
    )
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        accuracy_ax, action_ax = fig.axes[:2]

        legend = accuracy_ax.get_legend()
        assert legend is not None
        legend_box = legend.get_window_extent(renderer).padded(1)
        legend_collisions = [
            text.get_text()
            for text in accuracy_ax.texts
            if text.get_visible()
            and text.get_text().strip()
            and legend_box.overlaps(text.get_window_extent(renderer).padded(1))
        ]
        assert not legend_collisions, (
            f"panel-a legend overlaps annotations: {legend_collisions}"
        )

        error_only = next(
            line for line in accuracy_ax.lines if line.get_label() == "error-only RK4"
        )
        vertices = error_only.get_path().transformed(error_only.get_transform()).vertices
        tolerance_labels = [
            text for text in accuracy_ax.texts if r"\tau" in text.get_text()
        ]
        curve_collisions = []
        for text in tolerance_labels:
            text_box = text.get_window_extent(renderer).padded(1)
            intersects_curve = any(
                text_box.contains(
                    start[0] + (end[0] - start[0]) * fraction / 100,
                    start[1] + (end[1] - start[1]) * fraction / 100,
                )
                for start, end in zip(vertices, vertices[1:])
                for fraction in range(101)
            )
            if intersects_curve:
                curve_collisions.append(text.get_text())
        assert not curve_collisions, (
            f"panel-a tolerance labels intersect the error-only curve: {curve_collisions}"
        )

        bar_boxes = [
            patch.get_window_extent(renderer).padded(1)
            for container in action_ax.containers
            for patch in container.patches
            if patch.get_height() > 0
        ]
        bar_collisions = [
            text.get_text()
            for text in action_ax.texts
            if text.get_visible()
            and text.get_text().strip()
            and any(
                bar_box.overlaps(text.get_window_extent(renderer).padded(1))
                for bar_box in bar_boxes
            )
        ]
        assert not bar_collisions, (
            f"panel-b annotations overlap accepted-map bars: {bar_collisions}"
        )
    finally:
        plt.close(fig)


def test_component_regions_have_noncolor_encodings(
    monkeypatch, tmp_path: Path
) -> None:
    """Attached and detached regions must remain distinct in grayscale or color-deficient viewing."""
    generator = _load_generator()
    fig1 = _capture_generated_figure(
        generator, monkeypatch, tmp_path, "fig1_topology_conditioning"
    )
    fig3 = _capture_generated_figure(
        generator, monkeypatch, tmp_path, "fig2_robustness_noncommuting"
    )
    try:
        region_hatches = {
            collection.get_hatch()
            for collection in fig1.axes[0].collections
            if collection.get_label() in {"attached component", "detached component"}
        }
        assert len(region_hatches) == 2
        assert all(region_hatches)

        component_lines = [
            line
            for line in fig3.axes[0].lines
            if line.get_label() in {"attached component", "detached component"}
        ]
        assert len(component_lines) == 2
        assert len({line.get_linestyle() for line in component_lines}) == 2
        assert len({line.get_marker() for line in component_lines}) == 2
    finally:
        plt.close(fig1)
        plt.close(fig3)


def test_figure1_inset_has_an_explicit_zoom_cue(
    monkeypatch, tmp_path: Path
) -> None:
    """The reentrant-band inset must be explicit without obscuring panel-a evidence."""
    generator = _load_generator()
    fig = _capture_generated_figure(
        generator, monkeypatch, tmp_path, "fig1_topology_conditioning"
    )
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        main_ax = fig.axes[0]
        indicators = fig.findobj(match=InsetIndicator)
        assert len(indicators) == 1

        indicator = indicators[0]
        rectangle = indicator.rectangle
        assert rectangle.get_x() == pytest.approx(0.80)
        assert rectangle.get_y() == pytest.approx(0.0)
        assert rectangle.get_width() == pytest.approx(0.08)
        assert rectangle.get_height() == pytest.approx(2.9)

        visible_connectors = [
            connector
            for connector in indicator.connectors
            if connector is not None and connector.get_visible()
        ]
        assert len(visible_connectors) == 1

        connector = visible_connectors[0]
        connector_path = connector.get_path().transformed(
            connector.get_transform()
        )
        connector_vertices = connector_path.vertices
        connector_length = sum(
            float(((end - start) ** 2).sum() ** 0.5)
            for start, end in zip(
                connector_vertices[:-1], connector_vertices[1:]
            )
        )
        assert connector_length < 0.15 * main_ax.get_window_extent(renderer).width

        inset = next(
            axis for axis in main_ax.child_axes
            if axis.get_xlim() == pytest.approx((0.80, 0.88))
        )
        assert inset.patch.get_alpha() in (None, 1.0)
        assert inset.patch.get_facecolor()[:3] == pytest.approx((1.0, 1.0, 1.0))
        assert inset.spines["left"].get_linewidth() >= 0.9

        inset_box = inset.get_window_extent(renderer).padded(1)
        source_box = rectangle.get_window_extent(renderer).padded(1)
        source_tick_collisions = [
            label.get_text()
            for label in [*inset.get_xticklabels(), *inset.get_yticklabels()]
            if label.get_visible()
            and label.get_text().strip()
            and label.get_window_extent(renderer).overlaps(source_box)
        ]
        assert not source_tick_collisions, (
            "panel-a inset ticks overlap the highlighted source band: "
            f"{source_tick_collisions}"
        )

        hidden_components = []
        for collection in main_ax.collections:
            if collection.get_label() not in {
                "attached component",
                "detached component",
            }:
                continue
            intersects_component = any(
                path.transformed(collection.get_transform()).intersects_bbox(
                    inset_box, filled=True
                )
                for path in collection.get_paths()
            )
            if intersects_component:
                hidden_components.append(collection.get_label())
        assert not hidden_components, (
            f"panel-a inset hides component regions: {hidden_components}"
        )

        legend = main_ax.get_legend()
        assert legend is not None
        legend_box = legend.get_window_extent(renderer).padded(1)
        assert not inset_box.overlaps(legend_box)
        assert not connector_path.intersects_bbox(legend_box, filled=False)

        covered_annotations = [
            text.get_text()
            for text in main_ax.texts
            if text.get_visible()
            and text.get_text().strip()
            and inset_box.overlaps(text.get_window_extent(renderer).padded(1))
        ]
        assert not covered_annotations, (
            f"panel-a inset covers annotations: {covered_annotations}"
        )

        crossed_annotations = [
            text.get_text()
            for text in main_ax.texts
            if text.get_visible()
            and text.get_text().strip()
            and connector_path.intersects_bbox(
                text.get_window_extent(renderer).padded(1), filled=False
            )
        ]
        assert not crossed_annotations, (
            f"panel-a zoom leader crosses annotations: {crossed_annotations}"
        )

        assert inset.title.get_text() == ""
        inset_label = next(
            text for text in inset.texts if text.get_text() == "zoom"
        )
        label_box = inset_label.get_window_extent(renderer)
        assert inset_box.contains(*label_box.get_points()[0])
        assert inset_box.contains(*label_box.get_points()[1])
    finally:
        plt.close(fig)


def test_emitted_eps_uses_portable_embedded_glyphs(
    monkeypatch, tmp_path: Path
) -> None:
    """EPS text must survive Ghostscript conversion instead of disappearing."""
    generator = _load_generator()
    monkeypatch.setattr(generator, "FIG", tmp_path)

    with plt.rc_context(generator.JOURNAL_STYLE):
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "Arial EPS probe", ha="center")
        generator.save(fig, "eps_font_probe")

    eps = (tmp_path / "eps_font_probe.eps").read_text(
        encoding="latin-1"
    )
    assert "/FontType 3 def" in eps
    assert "/FontType 42 def" not in eps


def test_main_figures_do_not_depend_on_transparency(
    monkeypatch, tmp_path: Path
) -> None:
    """PDF and EPS variants must retain the same visible hierarchy."""
    generator = _load_generator()
    functions = (
        "fig1_topology_conditioning",
        "fig2_robustness_noncommuting",
        "fig3_candidate_maps",
        "fig4_certified_guard_conditioning",
        "fig5_tolerance_sweep",
    )

    for function_name in functions:
        fig = _capture_generated_figure(
            generator, monkeypatch, tmp_path, function_name
        )
        try:
            transparent = [
                (type(artist).__name__, artist.get_alpha())
                for artist in fig.findobj()
                if artist.get_visible()
                and artist.get_alpha() is not None
                and artist.get_alpha() < 1.0
            ]
            assert not transparent, (
                f"{function_name} relies on transparency that EPS cannot preserve: "
                f"{transparent[:8]}"
            )
        finally:
            plt.close(fig)
