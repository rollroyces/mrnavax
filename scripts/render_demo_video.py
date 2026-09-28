#!/usr/bin/env python3
"""Render the mrnavax demo output as a terminal-style MP4 video.

Reads the captured demo output from stdin (or a file) and renders each
line appearing line-by-line on a styled terminal canvas. Encodes to
MP4 via ffmpeg at 1280x720 @ 30 fps.

Usage:
        python examples/run_all.py | python scripts/render_demo_video.py > demo.mp4

Or, equivalently:
        ./render_demo_video.py examples/demo_capture.txt demo.mp4
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
# Style — matches the dark-theme mrnavax pipeline.svg (deep navy + cyan)
# ---------------------------------------------------------------------------
BG_COLOR = (16, 26, 46)  # #101a2e
FG_COLOR = (220, 230, 245)  # #dce6f5 — soft white for text
DIM_COLOR = (110, 125, 150)  # #6e7d96 — comments / metadata
ACCENT_COLOR = (93, 208, 255)  # #5dd0ff — labels + headers
GREEN_COLOR = (110, 220, 160)  # #6edca0 — tool names
PROMPT_COLOR = (155, 140, 255)  # #9b8cff — the "$" prompt
YELLOW_COLOR = (255, 200, 110)  # #ffc86e — "Done." success line

FRAME_W, FRAME_H = 1280, 720
FPS = 30

# Font: prefer the macOS SF Mono fallback chain.
FONT_CANDIDATES = [
    "/System/Library/Fonts/SFNSMono.ttf",
    "/System/Library/Fonts/Menlo.ttc",
    "/Library/Fonts/Menlo.ttc",
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf",
]


def _find_font(size: int) -> PIL.ImageFont.FreeTypeFont:
    for path in FONT_CANDIDATES:
        if Path(path).exists():
            return PIL.ImageFont.truetype(path, size)
    return PIL.ImageFont.load_default()


def _line_metrics(font: PIL.ImageFont.FreeTypeFont) -> tuple[int, int]:
    """Cell width + height in pixels for the given font."""
    # Use 'M' as a reference glyph for both width and full bounding box height.
    bbox = font.getbbox("M")
    cell_w = bbox[2] - bbox[0]
    ascent, descent = font.getmetrics()
    cell_h = ascent + descent
    return cell_w, cell_h


# ---------------------------------------------------------------------------
# Color picker — maps each demo line to a color based on its semantic role
# ---------------------------------------------------------------------------
def color_for(line: str) -> tuple[int, int, int]:
    """Pick a color for a line based on its semantic role in the demo output."""
    s = line.rstrip("\n")
    if not s.strip():
        return FG_COLOR
    if s.startswith("mrnavax "):
        return ACCENT_COLOR  # header
    if "Done." in s or "live AlphaGenome" in s or s.startswith("predict"):
        return YELLOW_COLOR  # success / highlight
    # Tool-output rows: '<tool_name>     <metric>'
    parts = s.split()
    if len(parts) >= 2 and parts[0] in {
        "codon", "neoantigen", "trial", "lnp", "scrna", "manufacture",
        "spatial", "construct", "utr-design", "variant-regulatory",
        "predict",
    }:
        return FG_COLOR
    if s.startswith("$") or s.startswith(">"):
        return PROMPT_COLOR
    if s.startswith("#") or s.startswith("//"):
        return DIM_COLOR
    return FG_COLOR


# ---------------------------------------------------------------------------
# Frame rendering
# ---------------------------------------------------------------------------
def render_frame(
    canvas: PIL.Image.Image,
    draw: PIL.ImageDraw.ImageDraw,
    font: PIL.ImageFont.FreeTypeFont,
    cell_w: int,
    cell_h: int,
    visible_lines: list[str],
    cursor_x: int,
    cursor_visible: bool,
    margin_x: int = 56,
    margin_y: int = 56,
) -> None:
    """Render one frame: terminal chrome + visible lines + blinking cursor."""
    # Background
    draw.rectangle([(0, 0), (FRAME_W, FRAME_H)], fill=BG_COLOR)

    # Title bar — "mrnavax · end-to-end demo"
    title_font = _find_font(18)
    draw.text(
        (margin_x, 18),
        "● mrnavax · end-to-end demo · ⌘",
        fill=DIM_COLOR,
        font=title_font,
    )

    # Separator line under title bar
    draw.line(
        [(margin_x, 48), (FRAME_W - margin_x, 48)],
        fill=(40, 50, 70),
        width=1,
    )

    # Visible lines
    for i, line in enumerate(visible_lines):
        y = margin_y + i * cell_h
        if y > FRAME_H - margin_y:
            break
        color = color_for(line)
        draw.text((margin_x, y), line, fill=color, font=font)

    # Blinking cursor
    if cursor_visible:
        cursor_y = margin_y + len(visible_lines) * cell_h
        if cursor_y <= FRAME_H - margin_y:
            draw.rectangle(
                [
                    (margin_x + cursor_x * cell_w, cursor_y),
                    (margin_x + cursor_x * cell_w + cell_w // 2, cursor_y + cell_h),
                ],
                fill=ACCENT_COLOR,
            )


def _encode_to_mp4(frame_dir: Path, output_path: Path, fps: int) -> None:
    """Pipe PNG frames to ffmpeg for H.264 encoding."""
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        raise RuntimeError("ffmpeg not found on PATH")

    pattern = str(frame_dir / "frame_%05d.png")
    cmd = [
        ffmpeg,
        "-y",
        "-framerate", str(fps),
        "-i", pattern,
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",  # ensure even dims
        "-preset", "medium",
        "-crf", "23",
        "-movflags", "+faststart",
        str(output_path),
    ]
    # NOTE: stderr goes to a file to avoid the 64KB deadlock
    # documented in the ascii-video references/troubleshooting.md
    log_path = frame_dir / "ffmpeg.log"
    with open(log_path, "w") as logf:
        result = subprocess.run(cmd, stderr=logf, stdout=subprocess.DEVNULL, check=False)
    if result.returncode != 0:
        msg = log_path.read_text() if log_path.exists() else "<no log>"
        raise RuntimeError(f"ffmpeg failed (rc={result.returncode}):\n{msg[-2000:]}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="Path to demo output text file, or '-' for stdin.",
    )
    parser.add_argument(
        "output",
        nargs="?",
        default="mrnavax_demo.mp4",
        help="Output MP4 path.",
    )
    parser.add_argument("--fps", type=int, default=FPS)
    parser.add_argument(
        "--hold-last-frames",
        type=int,
        default=45,
        help="Number of frames to hold the final state (default 45 = 1.5s).",
    )
    args = parser.parse_args()

    # Read input
    if args.input == "-":
        text = sys.stdin.read()
    else:
        text = Path(args.input).read_text()
    lines = text.splitlines()

    # Typewriter pacing: ~0.25s per line, hold for ~1.5s at the end
    FRAMES_PER_LINE = max(2, args.fps // 4)  # 0.25s per line

    font = _find_font(28)
    cell_w, cell_h = _line_metrics(font)
    print(f"[render_demo_video] font: cell_w={cell_w}, cell_h={cell_h}", file=sys.stderr)

    with tempfile.TemporaryDirectory(prefix="mrnavax_demo_") as tmp:
        frame_dir = Path(tmp)
        total_frames = len(lines) * FRAMES_PER_LINE + args.hold_last_frames
        print(
            f"[render_demo_video] {len(lines)} lines × {FRAMES_PER_LINE} frames/line "
            f"+ {args.hold_last_frames} hold = {total_frames} frames "
            f"({total_frames / args.fps:.1f}s @ {args.fps}fps)",
            file=sys.stderr,
        )

        # Pre-render to PNG
        for frame_idx in range(total_frames):
            canvas = PIL.Image.new("RGB", (FRAME_W, FRAME_H), BG_COLOR)
            draw = PIL.ImageDraw.Draw(canvas)

            line_idx = frame_idx // FRAMES_PER_LINE
            visible_lines = lines[: line_idx + 1]

            # Blinking cursor on the last visible line
            cursor_x = 0  # position 0 (column 0)
            # cursor blinks every ~15 frames (~0.5s)
            cursor_visible = (frame_idx % 30) < 15

            render_frame(
                canvas,
                draw,
                font,
                cell_w,
                cell_h,
                visible_lines,
                cursor_x,
                cursor_visible,
            )

            frame_path = frame_dir / f"frame_{frame_idx:05d}.png"
            canvas.save(frame_path, optimize=True)

        print(
            f"[render_demo_video] rendered {total_frames} PNG frames to {frame_dir}",
            file=sys.stderr,
        )

        # Encode to MP4
        _encode_to_mp4(frame_dir, Path(args.output), args.fps)
        print(f"[render_demo_video] wrote {args.output}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
