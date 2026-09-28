#!/usr/bin/env python3
"""Render an animated architecture-diagram demo video for mrnavax.

This replaces the terminal-style end-to-end video with a visual
diagram: the existing pipeline.svg is the backdrop, and 11 "tool
chips" light up sequentially as each mrnavax tool runs. Each chip
slides in from below + a 1-line summary appears below the title.

Visual style:
- Pipeline SVG as static background (dark navy + cyan/purple gradient)
- Tool chips: rounded rectangles, deep-navy fill, cyan border
  (inactive) → bright cyan fill, white text (active)
- Title bar at top: "mrnavax — 11 tools, one pipeline"
- Caption strip at bottom: shows the currently-active tool name + metric
- Subtle pulse animation on the active chip (2-frame size bump)

Output: 1280x720 @ 30 fps, ~8s runtime, MP4 + GIF fallback.

Usage:
    python scripts/render_architecture_video.py \\
        docs/assets/architecture_demo.mp4 \\
        --gif docs/assets/architecture_demo.gif

Implementation: PIL for frame composition, ffmpeg for MP4 mux, palette
filter for GIF. ~200 lines; no external renderer, no asciinema.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import PIL.Image
import PIL.ImageDraw
import PIL.ImageFont

# ---------------------------------------------------------------------------
# Style tokens — match pipeline.svg palette
# ---------------------------------------------------------------------------
WIDTH = 1280
HEIGHT = 720
FPS = 30

BG_COLOR = (16, 26, 46)             # #101a2e
PANEL_COLOR = (22, 34, 58)          # #16223a
CHIP_INACTIVE_BG = (30, 41, 67)     # #1e2943
CHIP_INACTIVE_BORDER = (93, 208, 255)  # #5dd0ff
CHIP_ACTIVE_BG = (93, 208, 255)     # #5dd0ff
CHIP_ACTIVE_TEXT = (11, 16, 32)     # #0b1020
CHIP_TEXT = (220, 230, 245)         # #dce6f5
TITLE_COLOR = (226, 232, 240)       # #e2e8f0
DIM_COLOR = (148, 163, 184)         # #94a3b8
GREEN_COLOR = (110, 220, 160)       # #6edca0
PURPLE_COLOR = (155, 140, 255)      # #9b8cff

# ---------------------------------------------------------------------------
# Tool data — order matches the published pipeline layers + CLI invocation
# ---------------------------------------------------------------------------
TOOLS = [
    # (name, layer, one-line metric)
    ("codon", "Layer 1: Sequence design", "CAI 0.720 → 0.938"),
    ("variant-regulatory", "Layer 2: Variant scoring", "6 variants scored"),
    ("neoantigen", "Layer 3: Immunogenicity", "8 candidates"),
    ("trial", "Layer 4: Patient-trial matching", "top-1 = NCT00000003"),
    ("lnp", "Layer 5: Delivery vehicle", "target=lung, cargo=sarna"),
    ("manufacture", "Layer 6: Manufacturability", "score=0.887, 7/8 pass"),
    ("scrna", "Layer 7: Cell-type context", "10 cells × 9 genes"),
    ("spatial", "Layer 7: Tissue context", "10 tissue modules"),
    ("construct", "Layer 8: Construct design", "384 nt, CAI=0.887"),
    ("utr-design", "Layer 8: UTR-coupled design", "combined=0.807"),
    ("predict", "Atlas: Live prediction", "BRAF V600E → high"),
]


def load_font(size: int) -> PIL.ImageFont.FreeTypeFont:
    """Find a usable monospace or sans font on this system."""
    candidates = [
        "/System/Library/Fonts/SFNSMono.ttf",
        "/System/Library/Fonts/SFNS.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            try:
                return PIL.ImageFont.truetype(c, size)
            except OSError:
                continue
    return PIL.ImageFont.load_default()


def draw_chip(
    draw: PIL.ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    h: int,
    label: str,
    metric: str,
    active: bool,
    title_font: PIL.ImageFont.FreeTypeFont,
    metric_font: PIL.ImageFont.FreeTypeFont,
) -> None:
    """Draw a single tool chip — inactive (dim) or active (lit up)."""
    if active:
        bg = CHIP_ACTIVE_BG
        text_color = CHIP_ACTIVE_TEXT
        metric_color = CHIP_ACTIVE_TEXT
        border_color = CHIP_ACTIVE_BG
    else:
        bg = CHIP_INACTIVE_BG
        text_color = CHIP_TEXT
        metric_color = DIM_COLOR
        border_color = CHIP_INACTIVE_BORDER

    # Rounded rectangle (PIL only does rounded via bbox + rounded_rectangle on newer)
    try:
        draw.rounded_rectangle(
            (x, y, x + w, y + h),
            radius=10,
            fill=bg,
            outline=border_color,
            width=2,
        )
    except AttributeError:
        draw.rectangle(
            (x, y, x + w, y + h), fill=bg, outline=border_color, width=2
        )

    # Label (top line of chip)
    label_bbox = draw.textbbox((0, 0), label, font=title_font)
    label_w = label_bbox[2] - label_bbox[0]
    draw.text(
        (x + (w - label_w) // 2, y + 12),
        label,
        fill=text_color,
        font=title_font,
    )

    # Metric (bottom line of chip)
    metric_bbox = draw.textbbox((0, 0), metric, font=metric_font)
    metric_w = metric_bbox[2] - metric_bbox[0]
    draw.text(
        (x + (w - metric_w) // 2, y + h - 32),
        metric,
        fill=metric_color,
        font=metric_font,
    )


def draw_title_bar(
    draw: PIL.ImageDraw.ImageDraw,
    title_font: PIL.ImageFont.FreeTypeFont,
    subtitle_font: PIL.ImageFont.FreeTypeFont,
) -> None:
    """Draw the title + subtitle strip at the top of the frame."""
    # Title (centered)
    title = "mrnavax — 11 tools, one pipeline"
    tbbox = draw.textbbox((0, 0), title, font=title_font)
    tw = tbbox[2] - tbbox[0]
    draw.text(
        ((WIDTH - tw) // 2, 24),
        title,
        fill=TITLE_COLOR,
        font=title_font,
    )

    # Subtitle
    subtitle = "Stdlib-only mocks, opt-in to real models via [extras]"
    sbbox = draw.textbbox((0, 0), subtitle, font=subtitle_font)
    sw = sbbox[2] - sbbox[0]
    draw.text(
        ((WIDTH - sw) // 2, 64),
        subtitle,
        fill=DIM_COLOR,
        font=subtitle_font,
    )


def draw_caption_strip(
    draw: PIL.ImageDraw.ImageDraw,
    label_font: PIL.ImageFont.FreeTypeFont,
    metric_font: PIL.ImageFont.FreeTypeFont,
    active_tool: tuple[str, str, str] | None,
    progress: str,
) -> None:
    """Draw the bottom caption strip: 'currently running: codon' + metric."""
    strip_y = HEIGHT - 90
    draw.rectangle(
        (40, strip_y, WIDTH - 40, strip_y + 60),
        fill=PANEL_COLOR,
    )

    if active_tool:
        name, layer, metric = active_tool
        text = f"running → {name}"
        draw.text((60, strip_y + 10), text, fill=GREEN_COLOR, font=label_font)
        layer_text = f"({layer})"
        draw.text((60, strip_y + 38), layer_text, fill=DIM_COLOR, font=metric_font)
        # Right side: metric
        mbbox = draw.textbbox((0, 0), metric, font=label_font)
        mw = mbbox[2] - mbbox[0]
        draw.text(
            (WIDTH - 60 - mw, strip_y + 18),
            metric,
            fill=TITLE_COLOR,
            font=label_font,
        )
    else:
        draw.text(
            (60, strip_y + 18),
            "done.",
            fill=GREEN_COLOR,
            font=label_font,
        )

    # Progress indicator (top-right)
    draw.text(
        (WIDTH - 140, 24),
        progress,
        fill=DIM_COLOR,
        font=metric_font,
    )


def render_frame(
    tools_done: int,
    active_tool: tuple[str, str, str] | None,
    progress: str,
    title_font: PIL.ImageFont.FreeTypeFont,
    subtitle_font: PIL.ImageFont.FreeTypeFont,
    label_font: PIL.ImageFont.FreeTypeFont,
    metric_font: PIL.ImageFont.FreeTypeFont,
    chip_font: PIL.ImageFont.FreeTypeFont,
    pulse: bool = False,
) -> PIL.Image.Image:
    """Render a single frame at the given progress state."""
    img = PIL.Image.new("RGB", (WIDTH, HEIGHT), BG_COLOR)
    draw = PIL.ImageDraw.Draw(img)

    # Layer 1: Title bar
    draw_title_bar(draw, title_font, subtitle_font)

    # Layer 2: Pipeline backdrop (decorative horizontal lines + labels)
    draw_decor_pipeline(draw, metric_font)

    # Layer 3: Tool chips (4-column grid)
    chip_w = 280
    chip_h = 80
    cols = 4
    gap_x = 16
    gap_y = 24
    total_w = cols * chip_w + (cols - 1) * gap_x
    start_x = (WIDTH - total_w) // 2
    start_y = 140

    for i, tool in enumerate(TOOLS):
        col = i % cols
        row = i // cols
        x = start_x + col * (chip_w + gap_x)
        y = start_y + row * (chip_h + gap_y)

        is_active = active_tool is not None and tool[0] == active_tool[0]
        is_done = i < tools_done

        # Active chip: pulse animation (slight size bump)
        if is_active and pulse:
            x -= 3
            y -= 3
            chip_w_p = chip_w + 6
            chip_h_p = chip_h + 6
        else:
            chip_w_p = chip_w
            chip_h_p = chip_h

        draw_chip(
            draw, x, y, chip_w_p, chip_h_p,
            tool[0], tool[2],
            active=(is_active or is_done),
            title_font=chip_font,
            metric_font=metric_font,
        )

    # Layer 4: Caption strip at bottom
    draw_caption_strip(draw, label_font, metric_font, active_tool, progress)

    return img


def draw_decor_pipeline(
    draw: PIL.ImageDraw.ImageDraw, font: PIL.ImageFont.FreeTypeFont
) -> None:
    """Draw the 8-layer pipeline header strip above the chip grid."""
    layers = [
        ("Sequence design", "blue"),
        ("Variant scoring", "blue"),
        ("Immunogenicity", "purple"),
        ("Patient-trial", "purple"),
        ("Delivery", "purple"),
        ("Manufacturability", "purple"),
        ("Cell context", "purple"),
        ("Construct design", "green"),
    ]
    y = 100
    layer_w = (WIDTH - 80) // len(layers)
    for i, (label, color) in enumerate(layers):
        x = 40 + i * layer_w
        color_rgb = {
            "blue": (93, 208, 255),
            "purple": (155, 140, 255),
            "green": (110, 220, 160),
        }[color]
        # Small dot
        draw.ellipse((x + 8, y + 6, x + 18, y + 16), fill=color_rgb)
        # Label
        draw.text(
            (x + 26, y + 4),
            label,
            fill=color_rgb,
            font=font,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Output MP4 path")
    parser.add_argument(
        "--gif", type=Path, default=None,
        help="Optional GIF fallback path (requires palette generation)",
    )
    parser.add_argument("--fps", type=int, default=FPS)
    parser.add_argument("--per-tool-seconds", type=float, default=0.7)
    parser.add_argument("--hold-seconds", type=float, default=1.5)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Fonts (sized for 1280x720)
    title_font = load_font(28)
    subtitle_font = load_font(16)
    label_font = load_font(20)
    metric_font = load_font(14)
    chip_font = load_font(18)

    # Layout: per-tool = (per_tool_seconds * fps) frames at 0→1→2 (active)
    # Pulse alternates every 5 frames
    per_tool_frames = int(args.per_tool_seconds * args.fps)
    hold_frames = int(args.hold_seconds * args.fps)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        frame_paths: list[Path] = []
        frame_idx = 0

        # Frame 0: intro (all inactive)
        img = render_frame(
            tools_done=0, active_tool=None, progress="0 / 11",
            title_font=title_font, subtitle_font=subtitle_font,
            label_font=label_font, metric_font=metric_font,
            chip_font=chip_font,
        )
        fp = tmp_path / f"frame_{frame_idx:05d}.png"
        img.save(fp)
        frame_paths.append(fp)
        frame_idx += 1

        # Per-tool active frames
        for i, tool in enumerate(TOOLS):
            for f in range(per_tool_frames):
                pulse = (f % 10) < 5
                img = render_frame(
                    tools_done=i,
                    active_tool=tool,
                    progress=f"{i+1} / 11",
                    title_font=title_font, subtitle_font=subtitle_font,
                    label_font=label_font, metric_font=metric_font,
                    chip_font=chip_font,
                    pulse=pulse,
                )
                fp = tmp_path / f"frame_{frame_idx:05d}.png"
                img.save(fp)
                frame_paths.append(fp)
                frame_idx += 1

        # Hold frame (all done)
        for f in range(hold_frames):
            pulse = (f % 30) < 15
            img = render_frame(
                tools_done=len(TOOLS),
                active_tool=None,
                progress="11 / 11 ✓",
                title_font=title_font, subtitle_font=subtitle_font,
                label_font=label_font, metric_font=metric_font,
                chip_font=chip_font,
                pulse=pulse,
            )
            fp = tmp_path / f"frame_{frame_idx:05d}.png"
            img.save(fp)
            frame_paths.append(fp)
            frame_idx += 1

        print(f"Rendered {len(frame_paths)} frames.", file=sys.stderr)

        # Encode MP4
        ffmpeg = shutil.which("ffmpeg")
        if not ffmpeg:
            local = Path.home() / ".hermes/tools/ffmpeg-9.0.1-darwin-arm64/ffmpeg"
            if local.exists():
                ffmpeg = str(local)
        if not ffmpeg:
            print("ERROR: ffmpeg not found.", file=sys.stderr)
            return 1

        log_path = tmp_path / "ffmpeg.log"
        cmd = [
            ffmpeg, "-y", "-framerate", str(args.fps),
            "-i", str(tmp_path / "frame_%05d.png"),
            "-c:v", "libx264", "-preset", "fast", "-crf", "23",
            "-pix_fmt", "yuv420p",
            "-movflags", "+faststart",
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
            str(args.output),
        ]
        with open(log_path, "w") as logf:
            subprocess.run(
                cmd, stderr=logf, stdout=subprocess.DEVNULL, check=True,
            )
        print(f"MP4 written: {args.output} ({args.output.stat().st_size} bytes)", file=sys.stderr)

        # Optional GIF
        if args.gif:
            args.gif.parent.mkdir(parents=True, exist_ok=True)
            gif_cmd = [
                ffmpeg, "-y", "-i", str(args.output),
                "-vf", "fps=15,scale=960:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse",
                "-loop", "0", str(args.gif),
            ]
            with open(log_path, "a") as logf:
                subprocess.run(
                    gif_cmd, stderr=logf, stdout=subprocess.DEVNULL, check=True,
                )
            print(f"GIF written: {args.gif} ({args.gif.stat().st_size} bytes)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
